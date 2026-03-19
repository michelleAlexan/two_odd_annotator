#%%

from ete4 import PhyloTree
from Bio import SeqIO
from pathlib import Path
import json
import random

from typing import Literal

from bio_tools.msa.mafft import run_mafft, trim_msa_by_gap_fraction
from bio_tools.phylo.fasttree import run_fasttree

from two_odd_annotator.utils.io import load_config
from two_odd_annotator.constants import (
    DEFAULT_CONFIG_PATH, 
    FILTERED_HMMER_FASTA, 
    FILTERED_DIAMOND_FASTA,
    FILTERED_BLASTP_FASTA, 
    ANNOTATION_FASTA,
    ANNOTATION_MSA, 
    ANNOTATION_MSA_TRIM,
    ANNOTATION_TREE
    )

# seed random for reproducibility
random.seed(42)

#load config
config = load_config(Path(__file__).parents[3] / DEFAULT_CONFIG_PATH)


# helper functions
def build_distance_lookup(tree):
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



#%%

# load test tree
test_tree_path = Path(__file__).parents[3] / "tests" / "data" / "test_anno_tree.nwk"
t = PhyloTree(open(test_tree_path), sp_naming_function=lambda name: name.split("__")[0].split("_")[-1])

dist_dict = build_distance_lookup(t)

# from the 354 leaves, randomly pick 25 leaves that will be considered candidates for annotation.
candidates = set()
# pick random indexes of the leaves list
random_indexes = random.sample(range(len(list(t.leaves()))), 25)

for idx in random_indexes:
    candidates.add(list(t.leaves())[idx].name)


#%%
def run(
        result_dir: Path, 
        config: dict, 
        seq_sim_method: Literal["hmmer", "diamond", "blastp", "all"] = "hmmer"
        ):

    major_minor_2ODD_clusters_path = config["annotate"]["major_minor_2ODD_clusters"]
    with open(major_minor_2ODD_clusters_path) as f:
        major_minor_2ODD_clusters = json.load(f)
    
    # merge ingroup 2ODDs and candidate genes into one fasta file
    ingroup_2ODD_fasta = config["annotate"]["ingroup"]
    if seq_sim_method == "hmmer":
        filtered_fasta = FILTERED_HMMER_FASTA
    elif seq_sim_method == "diamond":
        filtered_fasta = FILTERED_DIAMOND_FASTA
    elif seq_sim_method == "blastp":
        filtered_fasta = FILTERED_BLASTP_FASTA
    elif seq_sim_method == "all":
        # any file called filtered....fasta
        filtered_fasta = result_dir.glob("**/filtered_*.fasta")
    else:
        raise ValueError(f"Invalid seq_sim_method: {seq_sim_method}. Must be one of 'hmmer', 'diamond', 'blastp', or 'all'.")
    
    annotation_fasta_path = result_dir / ANNOTATION_FASTA

    seq_ids = set()
    with open(annotation_fasta_path, "w") as out_fasta:
        # write ingroup 2ODD sequences to annotation fasta
        for record in SeqIO.parse(ingroup_2ODD_fasta, "fasta"):
            if record.id in seq_ids:
                continue
            SeqIO.write(record, out_fasta, "fasta")
            seq_ids.add(record.id)

        # write candidate sequences that passed the pre-filtering step
        for subdir in result_dir.iterdir():
            if subdir.is_dir():
                for fasta_file in subdir.glob(filtered_fasta):
                    for record in SeqIO.parse(fasta_file, "fasta"):
                        if record.id in seq_ids:
                            continue
                        SeqIO.write(record, out_fasta, "fasta")
                        seq_ids.add(record.id)


    # run MSA
    annotation_msa_path = result_dir / ANNOTATION_MSA
    run_mafft(
        input_fasta=annotation_fasta_path, 
        output_fasta=annotation_msa_path
        )

    # trim MSA
    annotation_msa_trim_path = result_dir / ANNOTATION_MSA_TRIM
    trim_msa_by_gap_fraction(
        input_alignment=annotation_msa_path, 
        output_alignment=annotation_msa_trim_path, 
        gap_threshold=0.9
        )
    
    # run FastTree to build tree for annotation
    annotation_tree_path = result_dir / ANNOTATION_TREE
    run_fasttree(
        input_alignment=annotation_msa_trim_path, 
        output_tree=annotation_tree_path
        )
    
    # load tree
    tree = PhyloTree(str(annotation_tree_path), sp_naming_function=lambda name: name.split("__")[0].split("_")[-1])
    tax2names, tax2lineages, tax2rank = tree.annotate_ncbi_taxa(taxid_attr='species')
    tree.ladderize()  
    
     






# %%
