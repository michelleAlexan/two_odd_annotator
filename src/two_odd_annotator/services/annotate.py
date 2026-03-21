#%%

from ete4 import PhyloTree, Tree
from Bio import SeqIO
from pathlib import Path
import json

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


def get_candidate_headers(
        annotation_fasta: Path, 
        seq_to_2ODD_id: dict[str, str]
     ) -> set[str]:
    """
    Get all input sequences aka candidates that need to be annotated from the annotation fasta file.
    This is done by checking which sequence headers in the annotation fasta file are not present in the seq_to_2ODD_id dictionary,
    which contains the mapping of known 2ODD sequence headers to their corresponding 2ODD cluster IDs.
    """
    candidate_headers = set()
    for record in SeqIO.parse(annotation_fasta, "fasta"):
        if record.id not in seq_to_2ODD_id:
            candidate_headers.add(record.id)
    return candidate_headers


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
                      annotation_fasta: Path
                      ):
    """
    Take a ete4 Tree / PhyloTree and assign the 2ODD_id to leaf's properties. 
    If the leaf is an ingroup 2ODD, assign the corresponding ODD_id. If the leaf corresponds to a minor 2ODD cluster only, assig "minor_2ODD_cluster" as 2ODD_id.
    If the leaf is not an ingroup 2ODD, assign "candidate" as 2ODD_id.

    If the leaf is a characterized 2ODD, additionally assign "function" and "metabolic_pathway" to properties.

    Return modified tree. 
    """
    candidate_headers = get_candidate_headers(annotation_fasta, seq_to_2ODD_id)

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

    # Step 1: assign two_odd_ids without moving anything
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

            # assign
            if best_dist <= threshold and best_two_odd_id is not None:
                candidate.props["two_odd_id"] = best_two_odd_id
            else:
                candidate.props["two_odd_id"] = "unresolved"

    # Step 2: rebuild landscape in order
    resolved_landscape = []
    current_cluster = []
    current_id = None

    for cluster in landscape:
        for node in cluster:
            node_id = node.props["two_odd_id"]

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




#%%
def run(
        result_dir: Path, 
        config: dict, 
        seq_sim_method: Literal["hmmer", "diamond", "blastp", "all"] = "hmmer"
        ):
    
    annotation_fasta_path = result_dir / ANNOTATION_FASTA
    ingroup_2ODD_fasta = config["annotate"]["ingroup"]
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

    major_minor_2ODD_clusters_path = config["annotate"]["major_minor_2ODD_clusters"]
    with open(major_minor_2ODD_clusters_path) as f:
        major_minor_2ODD_clusters = json.load(f)






# %%
