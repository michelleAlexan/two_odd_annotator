#%%

from ete4 import PhyloTree, Tree
from Bio import SeqIO
from pathlib import Path
import json
import numpy as np
import pandas as pd

from typing import Literal, Callable

from bio_tools.msa.mafft import run_mafft, trim_msa_by_gap_fraction
from bio_tools.phylo.fasttree import run_fasttree
from bio_tools.viz.tree import (
    explore_tree_plant_groups,
    load_treecluster_assignments, 
    explore_tree_cluster_clades
    )


from two_odd_annotator.constants import (
    FILTERED_HMMER_FASTA, 
    FILTERED_DIAMOND_FASTA,
    FILTERED_BLASTP_FASTA, 
    ANNOTATION_FASTA,
    ANNOTATION_MSA, 
    ANNOTATION_MSA_TRIM,
    ANNOTATION_TREE
    )

#%% helper functions
def is_char_bait_sequence(leaf_name: str) -> bool:
    return len(leaf_name.split("__")) == 4

def reverse_major_minor_2ODD_dict(major_minor_2ODD_dict: dict[str, list[str]]) -> dict[str, str]:
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
    
    """
    seq_id_to_2ODD_id = {}

    for major_2ODD_id, seq_ids in major_minor_2ODD_dict.get("major_2ODDs", {}).items():
        for seq_id in seq_ids:
            seq_id_to_2ODD_id[seq_id] = major_2ODD_id

    for minor_2ODD_id, seq_ids in major_minor_2ODD_dict.get("minor_2ODDs", {}).items():
        for seq_id in seq_ids:
            seq_id_to_2ODD_id[seq_id] = minor_2ODD_id

    return seq_id_to_2ODD_id



def create_annotation_fasta(
    results_dir: Path,
    ingroup_2ODD_fasta: Path,
    output_fasta: Path,
    seq_sim_method: Literal["hmmer", "diamond", "blastp", "all"] = "hmmer"):

    # merge ingroup 2ODDs and candidate genes into one fasta file
    if seq_sim_method == "hmmer":
        filtered_fasta = FILTERED_HMMER_FASTA
    elif seq_sim_method == "diamond":
        filtered_fasta = FILTERED_DIAMOND_FASTA
    elif seq_sim_method == "blastp":
        filtered_fasta = FILTERED_BLASTP_FASTA
    elif seq_sim_method == "all":
        # any file called filtered....fasta
        filtered_fasta = results_dir.glob("**/filtered_*.fasta")
    else:
        raise ValueError(f"Invalid seq_sim_method: {seq_sim_method}. Must be one of 'hmmer', 'diamond', 'blastp', or 'all'.")
    

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
                for fasta_file in subdir.glob(filtered_fasta):
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
        tree_path: Path
        ):
    
    # run MSA
    run_mafft(
        input_fasta=fasta_path, 
        output_fasta=msa_path
        )
    
    # trim MSA
    trim_msa_by_gap_fraction(
        input_alignment=msa_path,
        output_alignment=msa_trim_path,
        gap_threshold=0.9
        )

    # run FastTree to build tree for annotation
    run_fasttree(
        input_alignment=msa_trim_path,
        output_tree=tree_path
        )
    
def assign_2ODD_props(tree: Tree, 
                      seq_to_2ODD_id: dict[str, str], 
                      candidate_headers: set[str]
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
            leaf.add_props(
                function=function,
                metabolic_pathway=metabolic_pathway
            )

        if leaf.name in seq_to_2ODD_id:
            two_odd_id = seq_to_2ODD_id[leaf.name]
            if "minor" in two_odd_id:
                leaf.add_props(two_odd_id="minor_2ODD_cluster")
            else:
                leaf.add_props(two_odd_id=two_odd_id)
        elif leaf.name in candidate_headers:
            leaf.add_props(two_odd_id="candidate")
        
        else:
            raise ValueError(f"Leaf name {leaf.name} not found in seq_to_2ODD_id or candidate_headers. \
                             Check if all sequences in the annotation fasta are accounted for in the seq_to_2ODD_id mapping or candidate_headers set.")


def assign_plant_group_props(tree:Tree |PhyloTree):
    """
    Take a ete4 Tree / PhyloTree and assign plant group properties to the leaves based on their taxonomic lineage.
        - Algae
        - Lycophytes (non-seed vascular plants)
        - Liverworts
        - Mosses
        - Ferns
        - Gymnosperms
        - Early Angiosperms (basal angiosperms)
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
            "Early Angiosperms": "#3BA0BC",
            "Monocots": "#7502d9",
            "Dicots": "#edc5ec",
        }
    
    def classify_plant(node):
        lineage = [t.lower() for t in node.props["named_lineage"]]

        if "zygnematophyceae" in lineage:
            return "Algae"
        elif "lycopodiopsida" in lineage:
            return "Lycophytes"    #Non-seed vascular plants
        elif "polypodiopsida" in lineage:
            return "Ferns"
        elif "bryophyta" in lineage or "anthocerotophyta" in lineage:
            return "Mosses"
        elif "marchantiophyta" in lineage:
            return "Liverworts"
        if "acrogymnospermae" in lineage:
            return "Gymnosperms"
        elif "liliopsida" in lineage:
            return "Monocots"
        elif any(x in lineage for x in ["eudicotyledons", 
                                        "magnoliopsida", 
                                        "mesangiospermae"]):
            return "Dicots"


        elif any(x in lineage for x in ["amborellales",
            "nymphaeales",
            "austrobaileyales", 
            "magnoliidae"]):
            return "Basal Angiosperms"
        else:
            print(f"Plant group couldnt be mapped for node {node.props['sci_name']}")
            print(lineage)

    for leaf in tree.leaves():
    
        plant_group = classify_plant(leaf)
        leaf.add_props(
            plant_group=plant_group,
            color=GROUP_COLORS.get(plant_group, None)
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
        names[i]: {
            names[j]: dm[i][j]
            for j in range(len(names))
        }
        for i in range(len(names))
    }

    return dist


def get_clusters(t: Tree | PhyloTree) -> list[list[Tree | PhyloTree]]:
    """
    Extract maximal monophyletic clusters based on 'two_odd_id',
    treating 'candidate' as transparent (ignored for identity, included in clusters).

    All leaves under a valid 2ODD clade are included, including edge candidates.
    """

    clusters = []
    node_to_id = {}

    # --- Step 1: bottom-up annotation (ignore candidates) ---
    for node in t.traverse("postorder"):
        if node.is_leaf:
            node_id = node.props.get("two_odd_id")
        else:
            child_ids = {
                node_to_id[child]
                for child in node.children
                if node_to_id[child] != "candidate"
            }

            if len(child_ids) == 1:
                node_id = next(iter(child_ids))  # pure real clade
            elif len(child_ids) == 0:
                node_id = "candidate"  # only candidates below
            else:
                node_id = None  # mixed real IDs

        node_to_id[node] = node_id

    assigned = set()

    # --- Step 2: extract maximal real clades ---
    for node in t.traverse("preorder"):
        node_id = node_to_id[node]

        # skip non-real clades
        if node_id is None or node_id == "candidate":
            continue

        parent = node.up
        parent_id = node_to_id.get(parent) if parent else None

        # maximal clade condition
        if node_id != parent_id:
            leaves = list(node.leaves())
            clusters.append(leaves)
            assigned.update(leaves)

    # --- Step 3: leftover candidates (not inside any real clade) ---
    for leaf in t.leaves():
        if leaf not in assigned:
            clusters.append([leaf])

    # --- Step 4: sort by leaf order ---
    leaf_order = {leaf: i for i, leaf in enumerate(t.leaves())}
    clusters.sort(key=lambda c: min(leaf_order[l] for l in c))

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
                dist_dict[a.name][b.name]
                for a in cluster_a
                for b in cluster_b
            )

            # strictly better OR tie → pick smaller index
            if (min_dist < best_dist) or (
                min_dist == best_dist and (best_j is None or j < best_j)
            ):
                best_dist = min_dist
                best_j = j

        neighbors[f"{i}"] = {
            "closest_cluster": best_j,
            "distance": round(best_dist, 3)
        }

    return neighbors



def two_odd_id_to_cluster_indices(clusters: list[list[Tree | PhyloTree]]) -> dict[str, list[int]]:
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
    (node.props.get("two_odd_id") for node in cluster if node.props.get("two_odd_id") != "candidate"),
    "candidates_only"
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
    candidates: set[str]
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
            (node.props.get("two_odd_id") for node in cluster if node.props.get("two_odd_id") != "candidate"),
            "candidates_only"
        )

        num_candidates = sum(1 for node in cluster if node.name in candidates)

        num_ingroup_2ODD = sum(
            1 for node in cluster
            if node.name in two_odd_ingroups and node.name not in candidates
        )

        # --- percentage ---
        percentage = None

        if (
            num_ingroup_2ODD > 0 and cluster_id in two_odd_to_ingroup_size
        ):
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




def clusters_meta_to_df(meta_info: dict[str, dict], neighboring_clusters: dict[int, tuple[int, float]]) -> pd.DataFrame:
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
            "neighboring_cluster_idx": int(neighboring_clusters[idx]["closest_cluster"]),
            "neighboring_cluster_dist": neighboring_clusters[idx]["distance"]
        }
        rows.append(row)

    df = pd.DataFrame(rows)

    # optional: sort by cluster index for readability
    df = df.sort_values("cluster_index").reset_index(drop=True)

    return df




#%%---------------------------------------------------------------
def run(
        result_dir: Path, 
        config: dict, 
        seq_sim_method: Literal["hmmer", "diamond", "blastp", "all"] = "hmmer"
        ):
    
    annotation_fasta_path = result_dir / ANNOTATION_FASTA
    ingroup_2ODD_fasta = config["annotate"]["ingroup"]
    major_minor_2ODD_dict = config["annotate"]["major_minor_2ODD_ids"]
    annotation_msa_path = result_dir / ANNOTATION_MSA
    annotation_msa_trim_path = result_dir / ANNOTATION_MSA_TRIM
    annotation_tree_path = result_dir / ANNOTATION_TREE


    create_annotation_fasta(
        results_dir=result_dir,
        ingroup_2ODD_fasta=ingroup_2ODD_fasta,
        output_fasta=annotation_fasta_path,
        seq_sim_method=seq_sim_method
    )


    from_fasta_to_nwk(
        fasta_path=annotation_fasta_path,
        msa_path=annotation_msa_path,
        msa_trim_path=annotation_msa_trim_path,
        tree_path=annotation_tree_path
    )

    # load tree
    tree = PhyloTree(str(annotation_tree_path), sp_naming_function=lambda name: name.split("__")[0].split("_")[-1])
    tax2names, tax2lineages, tax2rank = tree.annotate_ncbi_taxa(taxid_attr='species')
    tree.ladderize()  







# %%
