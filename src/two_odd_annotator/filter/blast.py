#%%
# todo: add logs folder in output directory to save all intermediate files and logs


import subprocess
import yaml
import pandas as pd
import os
from Bio import SeqIO

def load_config(path: str) -> dict:
    """Load YAML configuration file."""
    with open(path, "r") as yaml_file:
        config = yaml.safe_load(yaml_file)
    return config


def filtered_out_seqs(input_fasta, output_fasta):
    # Collect all sequence IDs from the input FASTA
    set_ids_input = {record.id for record in SeqIO.parse(input_fasta, "fasta")}
    # Collect all sequence IDs from the output FASTA
    set_ids_output = {record.id for record in SeqIO.parse(output_fasta, "fasta")}
    # Return IDs that are in input but not in output
    return set_ids_input.difference(set_ids_output)


def get_input_file_format():
    """Extract the format of a given input file and set a constant variable."""
    ... # todo

def get_num_sequences(fasta_file):
    """Count sequences in a FASTA file."""
    return sum(1 for _ in SeqIO.parse(fasta_file, "fasta"))


def run_blastp(input_file, output_dir, reference_db, evalue, threads, logs=True) -> pd.DataFrame:
    """
    Run BLASTP for one input file. 
    The blastp command takes an e-value, to return only those hits, 
    for which the found alignment score would occur maximum <evalue> amount of times
    just by chance in the given db. 

    Returns the Blast output file as pd.DataFrame (modified to include column names)
    """
    blast_output_file = os.path.join(output_dir, "blast_results.tsv")
    cmd = (
        f"blastp "
        f"-query {input_file} "
        f"-db {reference_db} "
        f"-out {blast_output_file} "
        f"-evalue {evalue} "
        f"-num_threads {threads} "
        f"-outfmt '6 qseqid sseqid pident length mismatch gapopen "
        f"qstart qend sstart send evalue bitscore'"
    )

    print(f"Run: {cmd}")
    result = subprocess.run(cmd, shell=True)
    if result.returncode != 0:
        raise RuntimeError(f"BLAST command failed: {cmd}")
    
    cols = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
    
    df = pd.read_csv(blast_output_file, sep="\t", header=None, names=cols)
 
    return df


def filter_by_thresholds(df, config, input_file, output_dir):
    """Filter BLAST results and write filtered FASTA file."""
    thresholds = config["tools"]["blast"]["thresholds"]
    pident = thresholds["pident"]
    length = thresholds["length"]
    bitscore = thresholds["bitscore"]
    num_hits = thresholds["num_hits"]

    df_filtered = (
        df[(df["pident"] >= pident) &
           (df["length"] >= length) &
           (df["bitscore"] >= bitscore)]
        .sort_values("bitscore", ascending=False)
        .groupby("qseqid")
        .head(num_hits)
        .reset_index(drop=True)
    )

    output_csv = os.path.join(output_dir, "blast_filtered_hits.csv")
    df_filtered.to_csv(output_csv, index=False)

    passed = set(df_filtered["qseqid"])
    output_fasta = os.path.join(output_dir, "filter_seq_sim_blast.fasta")

    count = 0
    with open(output_fasta, "w") as out:
        for record in SeqIO.parse(input_file, "fasta"):
            if record.id in passed:
                SeqIO.write(record, out, "fasta")
                count += 1

    print(f"  → {count} 2ODD sequences found (saved to {output_fasta})")
    return df_filtered



def main(input_file, config):
    """Run the full filter_seq_sim pipeline for a single FASTA file."""
    output_dir = config["output_dir"]

    # if enabled, check if the output has already been calculated and skip the computation
    if config["reuse_existing"]:
        blast_result_file = os.path.join(output_dir, "blast_results.tsv")
        filtered_hits = os.path.join(output_dir, "blast_filtered_hits.csv")
        filtered_fasta = os.path.join(output_dir, "filter_seq_sim_blast.fasta")
        if all(os.path.isfile(f) for f in [blast_result_file, filtered_hits, filtered_fasta]):
            print(f"SKIPPING {os.path.basename(output_dir)}: existing results found.")
            return 

    os.makedirs(output_dir, exist_ok=True)

    print(f"Processing {input_file}")
    print(f"  → input size: {get_num_sequences(input_file)} ")

    blast_cfg = config["tools"]["blast"]
    blast_output_df = run_blastp(
        input_file=input_file,
        output_dir=output_dir,
        reference_db=blast_cfg["reference_db"],
        evalue=blast_cfg["thresholds"]["evalue"],
        threads=blast_cfg["threads"],
    )


    filter_by_thresholds(blast_output_df, config, input_file, output_dir)
    print("\n")


def process_directory(input_dir, config):
    """Recursively find all FASTA files and run main() for each."""
    pep_fasta_files = [
        os.path.join(root, f)
        for root, _, files in os.walk(input_dir)
        for f in files 
        if f.lower().endswith((".fasta", ".fa", ".faa")) 
            and ".cds" not in f
    ]

    if not pep_fasta_files:
        raise ValueError(f"No FASTA files found in {input_dir}")

    print(f"Found {len(pep_fasta_files)} FASTA files.")
    counter = 0
    for f in pep_fasta_files:
        counter += 1
        print(f"----- process {counter} / {len(pep_fasta_files)} files-----")
        subdir = os.path.join(config["output_dir"], (os.path.basename(f).split(".")[0]))
        os.makedirs(subdir, exist_ok=True)
        # modify config for this file
        cfg_copy = config.copy()
        cfg_copy["output_dir"] = subdir
        main(f, cfg_copy)



#%%
if __name__ == "__main__":
    import time

    config = load_config("/Users/michellealexander/projects/FALGs/configs/filter.yml")
    input_path = config["input"]

    start_time = time.time()
    print(f"Starting filter_seq_sim - blast process at {time.ctime(start_time)}")

    if os.path.isfile(input_path):
        main(input_path, config)
    elif os.path.isdir(input_path):
        process_directory(input_path, config)
    else:
        raise ValueError("Input must be a FASTA file or a directory containing FASTA files.")

    end_time = time.time()
    print(f"Finished in {end_time - start_time:.1f} seconds")

# %%
