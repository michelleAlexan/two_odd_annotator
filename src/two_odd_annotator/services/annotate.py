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


def get_landscape(t: Tree | PhyloTree) -> list[list[Tree | PhyloTree]]:
    """
    Iterates over the leaves of the tree and groups them into sub-clusters based on their 'two_odd_id' property.

    Each tree node needs a 'two_odd_id' property.
    This can be a 2ODD cluster id, a minor 2ODD cluster id or cluster id is 'candidate'. 
    """
    leaves_list = list(t.leaves())

    landscape = []
    sub_cluster = [leaves_list[0]]

    # get two_odd_id of first leaf
    last_two_odd_id = leaves_list[0].props["two_odd_id"]


    for i in range(1, len(leaves_list)):
        leaf = leaves_list[i]
        current_two_odd_id = leaf.props["two_odd_id"]

        if current_two_odd_id == last_two_odd_id:
            sub_cluster.append(leaf)
        else:
            landscape.append(sub_cluster)
            sub_cluster = [leaf]
            last_two_odd_id = current_two_odd_id

    # append the last sub_cluster
    landscape.append(sub_cluster)


    return landscape

    
        
def resolve_candidates_in_landscape(
    landscape: list[list[Tree | PhyloTree]],
    dist_dict: dict[str, dict[str, float]],
    threshold: float
) -> list[list[Tree | PhyloTree]]:
    """
    Resolve candidate nodes in a linearized phylogenetic landscape (result of get_landscape function) 
    by assigning them to the closest neighboring 2ODD cluster, while preserving the original leaf order.


    Parameters
    ----------
    landscape : list[list[Tree | PhyloTree]]
        Ordered list of clusters. Each cluster is a list of tree nodes. Each node must have:
            - node.name : str
            - node.props["two_odd_id"] : str
        Cluster IDs are expected to be:
            - "candidate" which indicate user input sequences that need to be annotated
            - a 2ODD cluster ID (e.g. "2ODD01" or "mindor_2ODD_cluster")
    
    dist_dict : dict[str, dict[str, float]]
        Nested dictionary containing pairwise distances between sequences:
            dist_dict[seqA][seqB] = distance

    threshold : float
        Maximum distance allowed for assigning a candidate to a cluster.
        If the closest distance exceeds this threshold, the candidate is marked as "unresolved".

    Returns
    -------
    list[list[Tree | PhyloTree]]
        A new landscape where:
            - candidate nodes are assigned to neighboring clusters or marked as "unresolved"
            - clusters are rebuilt as contiguous blocks of nodes sharing the same two_odd_id
            - leaf order is preserved

    Algorithm
    ---------
    For each candidate node:
        1. Identify upstream cluster (previous in landscape, if any)
        2. Identify downstream cluster (next in landscape, if any)
        3. Compute distance to:
            - last node of upstream cluster
            - first node of downstream cluster
        4. Select the closest cluster
        5. If distance <= threshold:
               assign candidate to that cluster (via two_odd_id)
           else:
               assign "unresolved"

    After assignment:
        - Iterate through all nodes in original order
        - Group adjacent nodes with identical two_odd_id into new clusters

    Example
    -------
    Input landscape:
        [
            [candidate1, candidate2],
            [A1, A2],                  # two_odd_id = "2ODD01"
            [candidate3, candidate4],
            [B1, B2, B3]              # two_odd_id = "2ODD02"
        ]

    Suppose:
        - candidate2 is close to 2ODD01 → assigned to "2ODD01"
        - candidate3 is also closest to 2ODD01 → assigned to "2ODD01"
        - candidate4 is closest to 2ODD02 → assigned to "2ODD02"
        - candidate1 exceeds threshold → "unresolved"

    Resulting landscape (order preserved, clusters rebuilt):
        [
            [candidate1],                          # unresolved
            [candidate2, A1, A2, candidate3],      # 2ODD01
            [candidate4, B1, B2, B3]               # 2ODD02
        ]

    """

    # Step 1: assign targets without overwriting "candidate"
    for idx, cluster in enumerate(landscape):

        if cluster[0].props.get("two_odd_id") != "candidate":
            continue

        upstream_cluster = landscape[idx - 1] if idx > 0 else None
        downstream_cluster = landscape[idx + 1] if idx < len(landscape) - 1 else None

        for candidate in cluster:
            name = candidate.name

            best_two_odd_id = None
            best_dist = float("inf")

            # upstream
            if upstream_cluster:
                up_name = upstream_cluster[-1].name
                d = dist_dict[name][up_name]
                if d < best_dist:
                    best_dist = d
                    best_two_odd_id = upstream_cluster[0].props["two_odd_id"]

            # downstream
            if downstream_cluster:
                down_name = downstream_cluster[0].name
                d = dist_dict[name][down_name]
                if d < best_dist:
                    best_dist = d
                    best_two_odd_id = downstream_cluster[0].props["two_odd_id"]

            # assign WITHOUT overwriting candidate
            if best_dist <= threshold and best_two_odd_id is not None:
                candidate.props["_assigned_two_odd_id"] = best_two_odd_id
            else:
                candidate.props["two_odd_id"] = "unresolved"
                candidate.props["_assigned_two_odd_id"] = None


    # Step 2: rebuild landscape using assigned IDs (fallback to real ID)
    resolved_landscape = []
    current_cluster = []
    current_id = None

    for cluster in landscape:
        for node in cluster:
            node_id = (
                node.props.get("_assigned_two_odd_id")
                if node.props.get("two_odd_id") == "candidate"
                else node.props.get("two_odd_id")
            )

            if current_id is None:
                current_id = node_id
                current_cluster = [node]
                continue

            if node_id == current_id:
                current_cluster.append(node)
            else:
                resolved_landscape.append(current_cluster)
                current_cluster = [node]
                current_id = node_id

    if current_cluster:
        resolved_landscape.append(current_cluster)

    return resolved_landscape



def two_odd_id_to_landscape_indices(landscape: list[list[Tree | PhyloTree]]) -> dict[str, list[int]]:
    """
    Build a dictionary mapping each two_odd_id to the list of landscape indices where it occurs.

    Parameters
    ----------
    landscape : list[list[Tree | PhyloTree]]
        Ordered list of clusters. Each cluster is a list of tree nodes. Each node must have:
            - node.props["two_odd_id"] : str

    Returns
    -------
    dict[str, list[int]]
        Dictionary mapping each two_odd_id to the list of landscape indices where it occurs.
        e.g. {
            "2ODD01": [0, 1],
            "2ODD02": [2],
            "candidate": [3, 4, 5]
        }
    """
    id_to_indices = {}
    for idx, cluster in enumerate(landscape):
        two_odd_id = next(
    (node.props.get("two_odd_id") for node in cluster if node.props.get("two_odd_id") != "candidate"),
    None
)
        if two_odd_id not in id_to_indices:
            id_to_indices[two_odd_id] = []
        id_to_indices[two_odd_id].append(idx)
    return id_to_indices

def seq_id_to_landscape_idx(landscape: list[list[Tree | PhyloTree]]) -> dict[str, int]:
    """
    Build a dictionary mapping each sequence ID to its corresponding landscape index.

    Parameters
    ----------
    landscape : list[list[Tree | PhyloTree]]
        Ordered list of clusters. Each cluster is a list of tree nodes. Each node must have:
            - node.name : str

    Returns
    -------
    dict[str, int]
        Dictionary mapping each sequence ID to its corresponding landscape index.
        e.g. {
            "seqA": 0,
            "seqB": 1,
            "seqC": 2
        }
    """
    seq_id_to_idx = {}
    for idx, cluster in enumerate(landscape):
        for node in cluster:
            seq_id_to_idx[node.name] = idx
    return seq_id_to_idx



def landscape_meta_info(
    landscape: list[list[Tree | PhyloTree]],
    major_minor_2ODD_dict: dict[str, dict[str, list[str]]],
    two_odd_ingroups: set[str],
    candidates: set[str]
) -> dict[int, dict]:
    """
    Build metadata for each cluster in the landscape.
    """

    meta_info = {}

    # answering: how many sequences per 2ODD id were used as ingroup sequences for the annotation?
    two_odd_to_ingroup_size = {}
    for two_odd_id, seq_ids in major_minor_2ODD_dict.get("major_2ODDs", {}).items():
        two_odd_to_ingroup_size[two_odd_id] = len(
            {seq for seq in seq_ids if seq in two_odd_ingroups}
        )

    for idx, cluster in enumerate(landscape):

        cluster_id = next(
    (node.props.get("two_odd_id") for node in cluster if node.props.get("two_odd_id") != "candidate"),
    None
)
        size = len(cluster)

        is_unresolved = cluster_id == "unresolved"

        contains_candidates = any(node.name in candidates for node in cluster)

        # default
        percentage = None

        if (
            not is_unresolved
            and cluster_id in two_odd_to_ingroup_size
        ):
            
            # how many ingroup 2ODDs are in total in the tree?
            total = two_odd_to_ingroup_size[cluster_id]

            # how many ingroup sequences are in this landscape cluster?
            ingroup_count = sum(
                1 for node in cluster
                if node.name in two_odd_ingroups and node.name not in candidates
            )
            print(f"Cluster {idx} (2ODD id: {cluster_id}) has {ingroup_count} ingroup sequences out of {total} total ingroup sequences for that 2ODD cluster.")
            percentage = round(ingroup_count / total, 3)

        plant_groups = {
            node.props.get("plant_group")
            for node in cluster
            if node.props.get("plant_group") is not None
        }

        meta_info[idx] = {
            "two_odd_id": cluster_id,
            "percentage_of_ingroup_2ODD": percentage,
            "num_sequences": size,
            "plant_groups": sorted(plant_groups),
            "contains_candidates": contains_candidates,
        }

    return meta_info

def landscape_meta_to_df(meta_info: dict[int, dict]) -> pd.DataFrame:
    """
    Convert landscape metadata dictionary into a pandas DataFrame.

    Parameters
    ----------
    meta_info : dict[int, dict]
        Output from `landscape_meta_info`, where keys are cluster indices
        and values are metadata dictionaries.

    Returns
    -------
    pd.DataFrame
        DataFrame with one row per cluster and columns:
        - cluster_index
        - two_odd_id
        - percentage_of_ingroup_2ODD
        - num_sequences
        - plant_groups (comma-separated string)
        - contains_candidates
        - unresolved
    """

    rows = []

    for idx, data in meta_info.items():
        row = {
            "cluster_index": idx,
            "two_odd_id": data.get("two_odd_id"),
            "percentage_of_ingroup_2ODD": data.get("percentage_of_ingroup_2ODD"),
            "num_sequences": data.get("num_sequences"),
            "plant_groups": ", ".join(data.get("plant_groups", [])),
            "contains_candidates": data.get("contains_candidates"),
            "unresolved": data.get("unresolved"),
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
