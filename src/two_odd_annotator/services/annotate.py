# %% ======= IMPORTS =======

from ete4 import PhyloTree, Tree
from Bio import SeqIO
from pathlib import Path
import json
import pandas as pd
from collections import Counter

import re
from typing import Literal

from two_odd_annotator.utils.msa import run_mafft, trim_msa_by_gap_fraction
from two_odd_annotator.utils.phylo import run_fasttree


from two_odd_annotator.constants import (
    FILTERED_HMMER_FASTA,
    FILTERED_DIAMOND_FASTA,
    FILTERED_BLASTP_FASTA,
    ANNOTATION_FASTA,
    ANNOTATION_MSA,
    ANNOTATION_MSA_TRIM,
    ANNOTATION_TREE,
    ANNOTATION_CSV,
)

MAJOR_MINOR_2ODD_DICT = json.load(
    open(Path(__file__).parents[3] / "data/2ODDs/major_minor_2ODD_ids_manual.json")
)
MAJOR_2ODD_INFO_DICT = json.load(
    open(Path(__file__).parents[3] / "data/2ODDs/major_2ODD_char_info.json")
)


# %% ========= FUNCTIONS ========

def is_char_bait_sequence(leaf_name: str) -> bool:
    return len(leaf_name.split("__")) == 4


def reverse_major_minor_2ODD_dict(
    major_minor_2ODD_dict: dict[str, list[str]] = MAJOR_MINOR_2ODD_DICT,
) -> dict[str, str]:
    """
    Build a dictionary mapping each sequence ID to its corresponding 2ODD_id.

    Take the major_minor_2ODD_dict which has the format:
    {
        "major_2ODDs": {
            "2ODD01": [
                "seq_header1", ...],

        "minor_2ODDs": {
            "2ODD_minor14": [
                "seq_header2", ...],
             ...
        }
    }

    Return:
    {
        "seq_header1": "2ODD01",
        "seq_header2": "2ODD_minor14",
        ...
    }

    """
    seq_id_to_2ODD_id = {}

    for major_2ODD_id, seq_ids in major_minor_2ODD_dict.get("major_2ODDs", {}).items():
        for seq_id in seq_ids:
            seq_id_to_2ODD_id[seq_id] = major_2ODD_id

    for minor_2ODD_id, seq_ids in major_minor_2ODD_dict.get("minor_2ODDs", {}).items():
        for seq_id in seq_ids:
            seq_id_to_2ODD_id[seq_id] = minor_2ODD_id

    return seq_id_to_2ODD_id


SEQ_TO_2ODD_ID = reverse_major_minor_2ODD_dict(MAJOR_MINOR_2ODD_DICT)


def build_parent_map(two_odd_ids: list[str]) -> dict[str, str]:
    """
    Build parent map for nested 2ODD IDs.

    Assumes structure like:
        2ODD14
        2ODD14A
        2ODD14A1

    Returns:
        child -> immediate parent
    """

    id_set = set(two_odd_ids)
    parent_map = {}

    for tid in two_odd_ids:
        # extract base numeric part
        match = re.match(r"(2ODD\d+)", tid)
        if not match:
            continue

        base = match.group(1)

        if tid == base:
            continue  # top-level → no parent

        # try to find closest parent
        for i in range(len(tid) - 1, len(base) - 1, -1):
            candidate = tid[:i]
            if candidate in id_set:
                parent_map[tid] = candidate
                break

        # fallback: assign to base if present
        if tid not in parent_map and base in id_set:
            parent_map[tid] = base

    return parent_map


PARENT_DICT = build_parent_map(list(MAJOR_MINOR_2ODD_DICT["major_2ODDs"].keys()))


def get_base_id(two_odd_id: str | None):
    if two_odd_id is None or two_odd_id == "candidate":
        return two_odd_id
    return PARENT_DICT.get(two_odd_id, two_odd_id)


def create_annotation_fasta(
    results_dir: Path,
    ingroup_2ODD_fasta: Path,
    output_fasta: Path,
    seq_sim_method: Literal["hmmer", "diamond", "blastp", "all"] = "hmmer",
):

    # merge ingroup 2ODDs and candidate genes into one fasta file
    if seq_sim_method == "hmmer":
        filtered_fasta = FILTERED_HMMER_FASTA
    elif seq_sim_method == "diamond":
        filtered_fasta = FILTERED_DIAMOND_FASTA
    elif seq_sim_method == "blastp":
        filtered_fasta = FILTERED_BLASTP_FASTA
    elif seq_sim_method != "all":
        # allow the special "all" mode, which is handled below when
        # iterating over subdirectories; anything else is invalid
        raise ValueError(
            f"Invalid seq_sim_method: {seq_sim_method}. Must be one of 'hmmer', 'diamond', 'blastp', or 'all'."
        )

    # collect and write all pre-filtered sequences into one fasta file for annotation
    seq_ids = set()
    with open(output_fasta, "w") as out_fasta:
        # write ingroup 2ODD sequences to annotation fasta
        for record in SeqIO.parse(ingroup_2ODD_fasta, "fasta"):
            if record.id in seq_ids:
                continue
            SeqIO.write(record, out_fasta, "fasta")
            seq_ids.add(record.id)

        # write candidate sequences that passed the pre-filtering step
        for subdir in results_dir.iterdir():
            if subdir.is_dir():
                # when seq_sim_method == "all", include any filtered_*.fasta file;
                # otherwise, only include the file for the selected method
                if seq_sim_method == "all":
                    fasta_iter = subdir.glob("filtered_*.fasta")
                else:
                    fasta_iter = subdir.glob(filtered_fasta)

                for fasta_file in fasta_iter:
                    for record in SeqIO.parse(fasta_file, "fasta"):
                        if record.id in seq_ids:
                            continue
                        SeqIO.write(record, out_fasta, "fasta")
                        seq_ids.add(record.id)


def split_seqs_by_2ODD_membership(
    seq_to_2ODD_id: dict[str, str],
    tree: Tree | PhyloTree | None = None,
    fasta: Path | None = None,
) -> tuple[set[str], set[str]]:
    """
    Split sequence IDs into candidate vs ingroup based on presence in seq_to_2ODD_id.

    Supports either:
    - a phylogenetic tree (uses leaf names)
    - a FASTA file (uses record IDs)

    Parameters
    ----------
    seq_to_2ODD_id : dict[str, str]
        Mapping of known sequence IDs to 2ODD cluster IDs.

    tree : Tree | PhyloTree, optional
        Tree containing sequence leaves.

    fasta : Path, optional
        FASTA file containing sequences.

    Returns
    -------
    tuple[set[str], set[str]]
        (candidate_ids, ingroup_ids)

    Raises
    ------
    ValueError
        If neither or both of tree and fasta are provided.
    """

    if (tree is None and fasta is None) or (tree is not None and fasta is not None):
        raise ValueError("Provide exactly one of: tree OR fasta")

    candidate_ids = set()
    ingroup_ids = set()

    # Case 1: Tree
    if tree is not None:
        for leaf in tree.leaves():
            name = leaf.name
            if name in seq_to_2ODD_id:
                ingroup_ids.add(name)
            else:
                candidate_ids.add(name)

    # Case 2: FASTA
    else:
        for record in SeqIO.parse(fasta, "fasta"):
            if record.id in seq_to_2ODD_id:
                ingroup_ids.add(record.id)
            else:
                candidate_ids.add(record.id)

    return candidate_ids, ingroup_ids


def from_fasta_to_nwk(
    fasta_path: Path,
    msa_path: Path,
    msa_trim_path: Path,
    tree_path: Path,
    completed_msa: bool = False,
    completed_msa_trim: bool = False,
    completed_tree: bool = False,
):
    """
    From a FASTA file, run MSA and tree building to produce a Newick tree for annotation.
    If the MSA, trimmed MSA, or tree already exist (e.g. from a previous run), they will be reused.
    """

    # run MSA
    if not completed_msa:
        run_mafft(input_fasta=fasta_path, output_fasta=msa_path)

    # trim MSA
    if not completed_msa_trim:
        trim_msa_by_gap_fraction(
            input_alignment=msa_path, output_alignment=msa_trim_path, gap_threshold=0.9
        )

    # run FastTree to build tree for annotation
    if not completed_tree:
        run_fasttree(input_alignment=msa_trim_path, output_tree=tree_path)


def assign_2ODD_props(
    tree: Tree,
    seq_to_2ODD_id: dict[str, str] = SEQ_TO_2ODD_ID,
    candidate_headers: set[str] = set(),
):
    """
    Take a ete4 Tree / PhyloTree and assign the 2ODD_id to leaf's properties.
    If the leaf is an ingroup 2ODD, assign the corresponding ODD_id. If the leaf corresponds to a minor 2ODD cluster only, assig "minor_2ODD_cluster" as 2ODD_id.
    If the leaf is not an ingroup 2ODD, assign "candidate" as 2ODD_id.

    If the leaf is a characterized 2ODD, additionally assign "function" and "metabolic_pathway" to properties.

    Return modified tree.
    """

    for leaf in tree.leaves():

        if is_char_bait_sequence(leaf.name):
            accession, function, metabolic_pathway, tax_id = leaf.name.split("__")
            leaf.add_props(function=function, metabolic_pathway=metabolic_pathway)

        if leaf.name in seq_to_2ODD_id:
            two_odd_id = seq_to_2ODD_id[leaf.name]
            if "minor" in two_odd_id:
                leaf.add_props(two_odd_id="minor_2ODD_cluster")
            else:
                leaf.add_props(two_odd_id=two_odd_id)
        elif leaf.name in candidate_headers:
            leaf.add_props(two_odd_id="candidate")

        else:
            raise ValueError(
                f"Leaf name {leaf.name} not found in seq_to_2ODD_id or candidate_headers. \
                             Check if all sequences in the annotation fasta are accounted for in the seq_to_2ODD_id mapping or candidate_headers set."
            )


def assign_plant_group_props(tree: Tree | PhyloTree):
    """
    Take a ete4 Tree / PhyloTree and assign plant group properties to the leaves based on their taxonomic lineage.
        - Algae
        - Lycophytes (non-seed vascular plants)
        - Liverworts
        - Mosses
        - Ferns
        - Gymnosperms
        - Basal Angiosperms
        - Monocots
        - Dicots (eudicots)
    Note that the tree.annotate_ncbi_taxa() function must have been run beforehand to populate the "named_lineage" property for each leaf,
    which contains the taxonomic lineage as a list of taxonomic names.
    """

    GROUP_COLORS = {
        "Algae": "#574104",
        "Lycophytes": "#ab730c",
        "Liverworts": "#faaf00",
        "Ferns": "#c1d717",
        "Mosses": "#89be86",
        "Gymnosperms": "#076247",
        "Basal Angiosperms": "#3BA0BC",
        "Monocots": "#7502d9",
        "Dicots": "#edc5ec",
    }

    def classify_plant(node):
        lineage = [t.lower() for t in node.props["named_lineage"]]

        if "zygnematophyceae" in lineage:
            return "Algae"
        elif "lycopodiopsida" in lineage:
            return "Lycophytes"  # Non-seed vascular plants
        elif "polypodiopsida" in lineage:
            return "Ferns"
        elif "bryophyta" in lineage or "anthocerotophyta" in lineage:
            return "Mosses"
        elif "marchantiophyta" in lineage:
            return "Liverworts"
        if "acrogymnospermae" in lineage:
            return "Gymnosperms"
        elif any(
            x in lineage
            for x in ["amborellales", "nymphaeales", "austrobaileyales", "magnoliidae"]
        ):
            return "Basal Angiosperms"
        elif "liliopsida" in lineage:
            return "Monocots"
        elif any(
            x in lineage for x in ["eudicotyledons", "magnoliopsida", "mesangiospermae"]
        ):
            return "Dicots"

        else:
            print(f"Plant group couldnt be mapped for node {node.props['sci_name']}")
            print(lineage)

    for leaf in tree.leaves():
        plant_group = classify_plant(leaf)
        leaf.add_props(
            plant_group=plant_group, color=GROUP_COLORS.get(plant_group, None)
        )


def build_distance_lookup(tree):
    """
    Build a nested dictionary for pairwise distances between leaves in the tree.

    e.g. dist_dict[seqA][seqB] = distance
    """
    leaves = list(tree.leaves())
    names = [leaf.name for leaf in leaves]

    dm = tree.distance_matrix(squared=True)

    # build lookup dict
    dist = {
        names[i]: {names[j]: dm[i][j] for j in range(len(names))}
        for i in range(len(names))
    }

    return dist


def get_root_id(two_odd_id: str, parent_map: dict[str, str]) -> str:
    """
    Collapse nested 2ODD IDs to their root parent.
    Example: 11A -> 11, 11B -> 11
    """
    while two_odd_id in parent_map:
        two_odd_id = parent_map[two_odd_id]
    return two_odd_id


def get_clusters(
    t: Tree | PhyloTree, dist_dict: dict[str, dict[str, float]] = None,
    parent_map: dict[str, str] = PARENT_DICT,
) -> list[list[Tree | PhyloTree]]:
    """
    Extract maximal monophyletic clusters based on 'two_odd_id',
    treating 'candidate' as transparent (ignored for identity, included in clusters).

    All leaves under a valid 2ODD clade are included, including edge candidates.
    """

    def same_parent_only(child_ids: set[str]) -> bool:
        def root_id(x):
            return x.rstrip("ABCDEFGHIJKLMNOPQRSTUVWXYZ") if x else x

        roots = {root_id(cid) for cid in child_ids}
        return len(roots) == 1

    clusters = []
    node_to_id = {}

    # --- Step 1: bottom-up annotation (ignore candidates) ---
    for node in t.traverse("postorder"):
        if node.is_leaf:
            node_id = node.props.get("two_odd_id")
        else:
            child_ids = [
                node_to_id.get(child)
                for child in node.children
                if node_to_id.get(child) != "candidate"
            ]

            unique_ids = set(child_ids)

            if not unique_ids:
                node_id = "candidate"

            elif len(unique_ids) == 1:
                # PURE → keep original ID (including nested like 11B)
                node_id = next(iter(unique_ids))

            elif same_parent_only(unique_ids):
                # mixed but same parent → collapse to parent
                node_id = get_root_id(next(iter(unique_ids)), parent_map)

            else:
                node_id = None  # mixed families

        node_to_id[node] = node_id

    # --- Step 2: extract maximal real clades (top-down) ---
    assigned = set()

    for node in t.traverse("preorder"):
        node_id = node_to_id[node]

        if node_id is None:
            continue  # skip non-2ODD clades

        # skip nodes whose leaves are already fully assigned to a cluster
        if all(leaf in assigned for leaf in node.leaves()):
            continue

        parent = node.up
        parent_id = node_to_id.get(parent) if parent else None

        # check if this node is a valid cluster root (different from parent)
        # e.g., if node_id == 2ODD01 and parent_id == 2ODD01,
        #  then we are in the middle of a 2ODD01 clade and should not start a new cluster here.
        if node_id != parent_id:

            # all leaves under this node are part of the same cluster,
            # even if some are candidates or nested 2ODD IDs
            leaves = list(node.leaves())

            # identify nested and candidate leaves
            two_odd_id_leaves = []
            nested_2ODD_id_seq_to_id = {}
            candidate_leaves = []

            for leaf in leaves:
                leaf_2odd_id = node_to_id.get(leaf)
                if leaf_2odd_id == "candidate":
                    candidate_leaves.append(leaf)
                    continue

                # for example if leaf_2odd _id is 11A and node_id is 11,
                # then we consider this leaf as nested.
                # This is an edge case that needs to be handled
                # we know for example, that 2ODD14 represents F3H clade
                # and we also know that the FNSI clade sits inside the F3H clade.
                # Thus, 2ODD14A is a nested clade and we need to resolve this cluster separately.
                if (
                    leaf_2odd_id != node_id
                    and get_root_id(leaf_2odd_id, parent_map) == node_id
                ):
                    # the leaves under a 2ODD ID (e.g., 2ODD11) can contain multiple nested 2ODD IDs
                    # e.g. 2ODD11A and 2ODD11B
                    # this will be resolved below
                    nested_2ODD_id_seq_to_id[leaf] = leaf_2odd_id
                else:
                    # this list will never contain nested 2ODD IDs nor candidates,
                    # only "pure" 2ODD IDs that match the current node_id (e.g. 2ODD11)
                    two_odd_id_leaves.append(leaf)

            # --- there are no nested 2ODD IDs ---
            if not nested_2ODD_id_seq_to_id:
                # if no nested 2ODDs, simply assign all leaves (including candidates) to the same cluster
                clusters.append(leaves)
                assigned.update(leaves)
                continue

            # --- there are nested 2ODD IDs ---
            # makes things a bit more complicated
            # we want to separate the nested 2ODD IDs from the 2ODD ID above (ie., 2ODD14A from 2ODD14)
            # but we cannot just put them in a differnt cluster.
            # the nested cluster may be split into multiple clusters (e.g. not all share the same LCA)
            # and there can be candidate seqs that need to be resoved
            #
            # we must first identify the nested clusters
            # by finding the LCA of all leaves with the same nested 2ODD ID (e.g. 2ODD14A)
            # and then extract those nested clusters from the larger cluster (e.g. 2ODD14).
            assigned_nested_nodes = set()
            nested_clusters = []
            for subnode in node.traverse("preorder"):
                subnode_id = node_to_id.get(subnode)
                # still in the same cluster, e.g. 2ODD14, ignore
                if subnode_id == node_id:
                    continue  # same cluster, ignore

                if subnode in assigned_nested_nodes:
                    continue  # already assigned as part of another nested cluster, ignore

                # spotted a nested cluster, e.g. 2ODD14A
                elif get_root_id(subnode_id, parent_map) == node_id:
                    # this is a nested cluster (e.g. 2ODD14A inside 2ODD14)
                    assigned_nested_nodes.update(list(subnode.traverse("preorder")))

                    nested_leaves = list(
                        leaf
                        for leaf in subnode.leaves()
                        if leaf in nested_2ODD_id_seq_to_id
                    )
                    nested_clusters.append(nested_leaves)

            # --- resolve candidates ---
            # assign candidates to the closest nested cluster if distance dict is provided, otherwise assign to main cluster
            if dist_dict is not None and candidate_leaves:
                for candidate in candidate_leaves:
                    best_dist = float("inf")
                    best_nested_cluster_idx = None

                    for idx, nested_cluster in enumerate(nested_clusters):
                        for nested_leaf in nested_cluster:
                            dist = dist_dict[candidate.name][nested_leaf.name]
                            if dist < best_dist:
                                best_dist = dist
                                best_nested_cluster_idx = idx
                    for two_odd_id_leaf in two_odd_id_leaves:
                        dist = dist_dict[candidate.name][two_odd_id_leaf.name]
                        if dist < best_dist:
                            best_dist = dist
                            best_nested_cluster_idx = None  # assign to main cluster

                    # important: index 0 is a valid cluster index
                    if best_nested_cluster_idx is not None:
                        nested_clusters[best_nested_cluster_idx].append(candidate)
                    else:
                        two_odd_id_leaves.append(candidate)
            else:
                if candidate_leaves:
                    print(
                        f"Warning: candidate leaves found in cluster {node_id} but no distance dict provided."
                    )

            # --- store clusters ---
            if two_odd_id_leaves:
                clusters.append(two_odd_id_leaves)
                assigned.update(two_odd_id_leaves)
            if nested_clusters:
                clusters.extend(nested_clusters)
                assigned.update(leaf for cluster in nested_clusters for leaf in cluster)

    # --- Step 3: leftover candidates (not inside any real clade) ---
    for leaf in t.leaves():
        if leaf not in assigned:
            clusters.append([leaf])

    # --- Step 4: ensure leaves inside each cluster follow tree order,
    #             but keep the cluster list itself in traversal order
    leaf_order = {leaf: i for i, leaf in enumerate(t.leaves())}
    for cluster in clusters:
        cluster.sort(key=lambda l: leaf_order[l])

    return clusters


def compute_cluster_neighbors(clusters, dist_dict):
    neighbors = {}

    for i, cluster_a in enumerate(clusters):
        best_dist = float("inf")
        best_j = None

        for j, cluster_b in enumerate(clusters):
            if i == j:
                continue

            min_dist = min(
                dist_dict[a.name][b.name] for a in cluster_a for b in cluster_b
            )

            # strictly better OR tie → pick smaller index
            if (min_dist < best_dist) or (
                min_dist == best_dist and (best_j is None or j < best_j)
            ):
                best_dist = min_dist
                best_j = j

        neighbors[f"{i}"] = {"closest_cluster": best_j, "distance": round(best_dist, 3)}

    return neighbors


def two_odd_id_to_cluster_indices(
    clusters: list[list[Tree | PhyloTree]],
) -> dict[str, list[int]]:
    """
    Build a dictionary mapping each two_odd_id to the list of cluster indices where it occurs.

    Parameters
    ----------
    clusters : list[list[Tree | PhyloTree]]
        Ordered list of clusters. Each cluster is a list of tree nodes. Each node must have:
            - node.props["two_odd_id"] : str

    Returns
    -------
    dict[str, list[int]]
        Dictionary mapping each two_odd_id to the list of cluster indices where it occurs.
        e.g. {
            "2ODD01": [0, 1],
            "2ODD02": [2],
            "candidate": [3, 4, 5]
        }
    """
    id_to_indices = {}
    for idx, cluster in enumerate(clusters):
        two_odd_id = next(
            (
                node.props.get("two_odd_id")
                for node in cluster
                if node.props.get("two_odd_id") != "candidate"
            ),
            "candidates_only",
        )
        if two_odd_id not in id_to_indices:
            id_to_indices[two_odd_id] = []
        id_to_indices[two_odd_id].append(idx)
    return id_to_indices


def seq_to_cluster_idx(clusters: list[list[Tree | PhyloTree]]) -> dict[str, int]:
    """
    Build a dictionary mapping each sequence ID to its corresponding cluster index.

    Parameters
    ----------
    clusters : list[list[Tree | PhyloTree]]
        Ordered list of clusters. Each cluster is a list of tree nodes. Each node must have:
            - node.name : str

    Returns
    -------
    dict[str, int]
        Dictionary mapping each sequence ID to its corresponding cluster index.
        e.g. {
            "seqA": 0,
            "seqB": 1,
            "seqC": 2
        }
    """
    seq_id_to_idx = {}
    for idx, cluster in enumerate(clusters):
        for node in cluster:
            seq_id_to_idx[node.name] = idx
    return seq_id_to_idx


def cluster_meta_info(
    clusters: list[list[Tree | PhyloTree]],
    major_minor_2ODD_dict: dict[str, dict[str, list[str]]],
    two_odd_ingroups: set[str],
    candidates: set[str],
) -> dict[int, dict]:
    """
    Build metadata for each cluster in the clusters.
    """

    meta_info = {}

    # total ingroup size per 2ODD
    two_odd_to_ingroup_size = {}
    for two_odd_id, seq_ids in major_minor_2ODD_dict.get("major_2ODDs", {}).items():
        two_odd_to_ingroup_size[two_odd_id] = len(
            {seq for seq in seq_ids if seq in two_odd_ingroups}
        )

    for idx, cluster in enumerate(clusters):

        cluster_id = next(
            (
                node.props.get("two_odd_id")
                for node in cluster
                if node.props.get("two_odd_id") != "candidate"
            ),
            "candidates_only",
        )

        num_candidates = sum(1 for node in cluster if node.name in candidates)

        num_ingroup_2ODD = sum(
            1
            for node in cluster
            if node.name in two_odd_ingroups and node.name not in candidates
        )

        # --- percentage ---
        percentage = None

        if num_ingroup_2ODD > 0 and cluster_id in two_odd_to_ingroup_size:
            total = two_odd_to_ingroup_size[cluster_id]

            percentage = round(num_ingroup_2ODD / total, 3)

        # --- plant groups ---
        plant_groups = {
            node.props.get("plant_group")
            for node in cluster
            if node.props.get("plant_group") is not None
        }

        meta_info[f"{idx}"] = {
            "two_odd_id": cluster_id,
            "percentage_of_ingroup_2ODD": percentage,
            "num_ingroup_2ODD": num_ingroup_2ODD,
            "num_candidates": num_candidates,
            "plant_groups": sorted(plant_groups),
        }

    return meta_info


def clusters_meta_to_df(
    meta_info: dict[str, dict], neighboring_clusters: dict[int, tuple[int, float]]
) -> pd.DataFrame:
    """
    Convert clusters metadata dictionary into a pandas DataFrame.

    Parameters
    ----------
    meta_info : dict[str, dict]
        Output from `clusters_meta_info`, where keys are cluster indices
        and values are metadata dictionaries.

    neighboring_clusters : dict[int, tuple[int, float]]
        Dictionary mapping each cluster index to a tuple of its neighboring cluster index and distance.

    Returns
    -------
    pd.DataFrame
        DataFrame with one row per cluster and columns:
        - cluster_index
        - two_odd_id
        - perc_of_ingroup_2ODD
        - n_ingroup_2ODD
        - n_candidates
        - plant_groups (comma-separated string)
        - neighboring_cluster_idx
        - neighboring_cluster_dist
    """

    rows = []

    for idx, data in meta_info.items():
        row = {
            "cluster_index": int(idx),
            "two_odd_id": data.get("two_odd_id"),
            "perc_of_ingroup_2ODD": data.get("percentage_of_ingroup_2ODD"),
            "n_ingroup_2ODD": data.get("num_ingroup_2ODD"),
            "n_candidates": data.get("num_candidates"),
            "plant_groups": ", ".join(data.get("plant_groups", [])),
            "neighboring_cluster_idx": int(
                neighboring_clusters[idx]["closest_cluster"]
            ),
            "neighboring_cluster_dist": neighboring_clusters[idx]["distance"],
        }
        rows.append(row)

    df = pd.DataFrame(rows)

    # optional: sort by cluster index for readability
    df = df.sort_values("cluster_index").reset_index(drop=True)

    return df


def get_candidate_to_char_baits_df(
    candidate_headers, ingroup_headers, dist_dict, seq_id_to_idx_dict
):
    """
    Assign each candidate sequence to its closest characterized bait sequences.

    Parameters
    ----------
    candidate_headers : list of str
        Candidate sequence headers.
    ingroup_headers : list of str
        All sequence headers in the ingroup tree.
    dist_dict : dict of dict
        dist_dict[candidate][bait] gives the distance between candidate and bait.
    seq_id_to_idx_dict : dict
        Mapping of sequence ID to cluster index.

    Returns
    -------
    pd.DataFrame
        Each row is a candidate sequence with cluster index, closest and second closest
        characterized bait, and corresponding distances.
    """
    # get characterized baits
    char_baits = [h for h in ingroup_headers if is_char_bait_sequence(h)]

    candidate_char_baits_dict = {}
    for candidate in candidate_headers:
        # sort closest baits by distance
        closest_char_baits = sorted(
            char_baits, key=lambda bait: dist_dict[candidate][bait]
        )[:3]

        candidate_char_baits_dict[candidate] = {
            "cluster_idx": seq_id_to_idx_dict.get(candidate, None),
            "closest_char_bait": closest_char_baits[0] if closest_char_baits else None,
            "closest_char_bait_dist": (
                dist_dict[candidate][closest_char_baits[0]]
                if closest_char_baits
                else None
            ),
            "second_closest_char_bait": (
                closest_char_baits[1] if len(closest_char_baits) > 1 else None
            ),
            "second_closest_char_bait_dist": (
                dist_dict[candidate][closest_char_baits[1]]
                if len(closest_char_baits) > 1
                else None
            ),
        }

    # If there are no candidates, return an empty DataFrame with the
    # expected schema so downstream merges on "cluster_idx" work.
    if not candidate_char_baits_dict:
        return pd.DataFrame(
            columns=[
                "candidate",
                "cluster_idx",
                "closest_char_bait",
                "closest_char_bait_dist",
                "second_closest_char_bait",
                "second_closest_char_bait_dist",
            ]
        )

    df = (
        pd.DataFrame.from_dict(candidate_char_baits_dict, orient="index")
        .reset_index()
        .rename(columns={"index": "candidate"})
    )
    return df


def add_consensus_annotations(df, threshold=0.9, min_size=2):
    """
    Add consensus function and pathway columns based on characterized bait sequences.

    Parameters
    ----------
    df : pd.DataFrame
        Must contain column 'associated_characterized_bait_sequences' (list of strings).
    threshold : float
        Fraction required for consensus (default: 0.9).
    min_size : int
        Minimum number of sequences required (default: 2).

    Returns
    -------
    pd.DataFrame
        With added columns:
        - consensus_function
        - consensus_pathway
    """

    def get_consensus(values):
        if len(values) < min_size:
            return None

        counts = Counter(values)
        most_common, count = counts.most_common(1)[0]

        if count / len(values) >= threshold:
            return most_common
        return None

    consensus_functions = []
    consensus_pathways = []

    for ids in df["associated_characterized_bait_sequences"]:
        if not ids:
            consensus_functions.append(None)
            consensus_pathways.append(None)
            continue

        # extract function + pathway
        functions = [seq.split("__")[1] for seq in ids if "__" in seq]
        pathways = [seq.split("__")[2] for seq in ids if "__" in seq]

        consensus_functions.append(get_consensus(functions))
        consensus_pathways.append(get_consensus(pathways))

    df = df.copy()
    df["consensus_function"] = consensus_functions
    df["consensus_pathway"] = consensus_pathways

    return df


def add_annotation_columns(df, threshold=0.8):
    df = df.copy()

    # --- condition: cluster is well resolved ---
    well_resolved = df["perc_of_ingroup_2ODD"] >= threshold

    # --- annotations ---
    df["annotated_two_odd_id"] = df["two_odd_id"].where(well_resolved, None)

    df["annotated_function"] = df["consensus_function"].where(
        well_resolved & df["consensus_function"].notna(), None
    )

    df["annotated_metabolic_pathway"] = df["consensus_pathway"].where(
        well_resolved & df["consensus_pathway"].notna(), None
    )

    # --- REORDER COLUMNS ---
    first_cols = [
        "candidate",
        "annotated_two_odd_id",
        "annotated_function",
        "annotated_metabolic_pathway",
    ]

    cluster_cols = [
        "cluster_index",
        "two_odd_id",
        "perc_of_ingroup_2ODD",
        "n_ingroup_2ODD",
        "n_candidates",
        "plant_groups",
        "neighboring_cluster_idx",
        "neighboring_cluster_dist",
    ]

    # everything else = bait + metadata columns
    other_cols = [col for col in df.columns if col not in first_cols + cluster_cols]

    ordered_cols = first_cols + cluster_cols + other_cols

    return df[ordered_cols]


# %%---------------------------------------------------------------
def run(
    result_dir: Path,
    config: dict,
    seq_sim_method: Literal["hmmer", "diamond", "blastp", "all"] = "hmmer",
    completed_annotation_steps: dict | None = None,
):

    annotation_fasta_path = result_dir / ANNOTATION_FASTA
    ingroup_2ODD_fasta = config["annotate"]["bait_sequence_collection"]
    major_minor_2ODD_dict = json.load(open(config["annotate"]["major_minor_2ODD_ids"]))
    major_char_2ODD_info_dict = json.load(
        open(config["annotate"]["major_2ODDs_functional_characterization"])
    )
    annotation_msa_path = result_dir / ANNOTATION_MSA
    annotation_msa_trim_path = result_dir / ANNOTATION_MSA_TRIM
    annotation_tree_path = result_dir / ANNOTATION_TREE
    save_tree = config["annotate"].get("save_tree", False)
    reuse_existing = config["pipeline"].get("reuse_existing", False)

    if not reuse_existing:
        completed_annotation_steps = (
            None  # ignore already completed steps if not reusing existing results
        )

    # ========== FROM FASTA TO NWK ==========

    if (
        completed_annotation_steps is None
        or not completed_annotation_steps["annotation_fasta"]
    ):
        create_annotation_fasta(
            results_dir=result_dir,
            ingroup_2ODD_fasta=ingroup_2ODD_fasta,
            output_fasta=annotation_fasta_path,
            seq_sim_method=seq_sim_method,
        )

    if (
        completed_annotation_steps is None
        or not completed_annotation_steps["annotation_tree"]
    ):
        completed_msa = (
            completed_annotation_steps.get("annotation_msa")
            if completed_annotation_steps
            else False
        )
        completed_msa_trim = (
            completed_annotation_steps.get("annotation_msa_trim")
            if completed_annotation_steps
            else False
        )
        completed_tree = (
            completed_annotation_steps.get("annotation_tree")
            if completed_annotation_steps
            else False
        )

        from_fasta_to_nwk(
            fasta_path=annotation_fasta_path,
            msa_path=annotation_msa_path,
            msa_trim_path=annotation_msa_trim_path,
            tree_path=annotation_tree_path,
            completed_msa=completed_msa,
            completed_msa_trim=completed_msa_trim,
            completed_tree=completed_tree,
        )

    # load tree
    tree = PhyloTree(
        open(annotation_tree_path), sp_naming_function=lambda name: name.split("__")[-1]
    )
    tax2names, tax2lineages, tax2rank = tree.annotate_ncbi_taxa(taxid_attr="species")
    tree.ladderize()
    if save_tree:
        tree.write(outfile=str(result_dir / ANNOTATION_TREE))

    # ========== ANNOTATE TREE WITH 2ODD IDS AND PLANT GROUPS ==========
    ingroup_seq_to_2ODD_id = reverse_major_minor_2ODD_dict(major_minor_2ODD_dict)
    candidate_headers, ingroup_headers = split_seqs_by_2ODD_membership(
        seq_to_2ODD_id=ingroup_seq_to_2ODD_id, tree=tree
    )
    assign_plant_group_props(tree)

    assign_2ODD_props(
        tree=tree,
        seq_to_2ODD_id=ingroup_seq_to_2ODD_id,
        candidate_headers=candidate_headers,
    )

    # ========== CREATE DATAFRAMES ==========

    # cluster_df:
    # representing the clusters of sequences in the tree,
    # their assigned 2ODD ID, the percentage of ingroup 2ODDs in the cluster,
    # and the plant groups represented in the cluster
    dist_dict = build_distance_lookup(tree)
    clusters = get_clusters(tree, dist_dict=dist_dict)
    neighboring_clusters = compute_cluster_neighbors(clusters, dist_dict)
    meta_info = cluster_meta_info(
        clusters=clusters,
        major_minor_2ODD_dict=major_minor_2ODD_dict,
        two_odd_ingroups=ingroup_headers,
        candidates=candidate_headers,
    )
    cluster_df = clusters_meta_to_df(meta_info, neighboring_clusters)

    # 2ODD_info_df:
    # all 2ODD ids and, if applicable, their associated functions and metabolic pathways shown in experiments
    two_odd_info_df = (
        pd.DataFrame.from_dict(major_char_2ODD_info_dict, orient="index")
        .reset_index()
        .rename(columns={"index": "two_odd_id"})
    )
    two_odd_info_df.rename(
        columns={
            "functions": "associated_functions",
            "metabolic_pathways": "associated_metabolic_pathways",
            "char_2ODD_ids": "associated_characterized_bait_sequences",
        },
        inplace=True,
    )
    two_odd_info_df = add_consensus_annotations(two_odd_info_df)

    # candidate_char_baits_df:
    # each row is a candidate sequence
    # and gets assigned cluster idx, closest, and second closest characterized 2ODD sequence
    # and the corresponding distance
    seq_id_to_idx = seq_to_cluster_idx(clusters)

    candidate_char_baits_df = get_candidate_to_char_baits_df(
        candidate_headers=candidate_headers,
        ingroup_headers=ingroup_headers,
        dist_dict=dist_dict,
        seq_id_to_idx_dict=seq_id_to_idx,
    )

    # create left join of candidate_char_baits_df with two_odd_info_df and cluster_df to get a full picture of each candidate's assigned cluster, closest characterized baits, and the known functions of those baits
    candidate_char_baits_full_df = candidate_char_baits_df.merge(
        cluster_df, how="left", left_on="cluster_idx", right_on="cluster_index"
    ).merge(two_odd_info_df, how="left", left_on="two_odd_id", right_on="two_odd_id")
    annotation_df = add_annotation_columns(candidate_char_baits_full_df)

    # save dataframe
    annotation_df.to_csv(result_dir / ANNOTATION_CSV, index=False)

    return annotation_df


# %%
