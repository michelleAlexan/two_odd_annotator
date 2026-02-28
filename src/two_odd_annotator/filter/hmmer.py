#%%
import subprocess
import yaml
import pandas as pd
import os
from Bio import SeqIO
from pathlib import Path
import re

#%%

def get_accessions_safe(df, column="Sequence"):
    """
    Safely extract a set of accessions from a dataframe column.
    
    Handles:
      - df is None
      - df has zero rows
      - df missing the column
      - df loaded without headers

    Returns
    -------
    set
        A set of accessions (may be empty).
    """

    # Case 1 — df is None
    if df is None:
        return set()

    # Case 2 — not a DataFrame (e.g., user passed wrong type)
    if not isinstance(df, pd.DataFrame):
        return set()

    # Case 3 — empty dataframe (no rows)
    if df.empty:
        return set()

    # Case 4 — missing expected column
    if column not in df.columns:
        return set()

    # Normal case
    return set(df[column].dropna().astype(str))


def filter_fasta_by_headers(input_file, output_dir, accession_list):
    """
    Write a new FASTA file containing only sequences whose FASTA headers
    (the part after '>' up to the first space) exactly match items in accession_list.
    """
    input_fasta = Path(input_file)
    output_fasta = Path(os.path.join(output_dir, "filter_seq_sim_hmmer.fasta")) 
    keep = set(accession_list)

    with input_fasta.open() as fin, output_fasta.open("w") as fout:
        write_seq = False
        for line in fin:
            if line.startswith(">"):
                header = line[1:].strip().split()[0]   # extract accession
                write_seq = header in keep
                if write_seq:
                    fout.write(line)
            else:
                if write_seq:
                    fout.write(line)


def parse_hmmsearch_output(hmmout_path):
    """
    Parse the 'Scores for complete sequences' table from an hmmsearch output file
    and return a pandas DataFrame.

    Parameters
    ----------
    hmmout_path : str or Path
        Path to an hmmsearch output text file.

    Returns
    -------
    pandas.DataFrame
        Columns:
        ['full_Evalue','full_score','full_bias',
         'bestdom_Evalue','bestdom_score','bestdom_bias',
         'exp','N','Sequence','Description']
    """

    rows = []
    in_table = False

    with open(hmmout_path) as f:
        for line in f:
            # Identify the table header that signals the start
            if "Scores for complete sequences" in line:
                in_table = False  # header appears before table
            if re.match(r"\s+E-value\s+score", line):
                in_table = True   # the real table starts after this
                continue
            if "---" in line and in_table:  # skip the dashed header line
                continue

            # Stop parsing if table ends (blank line after rows)
            if in_table and line.strip() == "":
                break

            # Parse table rows
            if in_table:
                # Split into columns **up to Sequence**, then keep the rest as Description
                parts = line.rstrip("\n").split()
                if len(parts) < 9:
                    continue  # skip malformed rows

                (
                    full_eval, full_score, full_bias,
                    best_eval, best_score, best_bias,
                    exp, N, seq_id,
                ) = parts[:9]

                # Everything after the first 9 fields = description
                description = " ".join(parts[9:]) if len(parts) > 9 else ""

                rows.append(
                    {
                        "full_Evalue": full_eval,
                        "full_score": float(full_score),
                        "full_bias": float(full_bias),
                        "bestdom_Evalue": best_eval,
                        "bestdom_score": float(best_score),
                        "bestdom_bias": float(best_bias),
                        "exp": float(exp),
                        "N": int(N),
                        "Sequence": seq_id,
                        "Description": description,
                    }
                )

    return pd.DataFrame(rows)

def filter_hmmer_df(df:pd.DataFrame) -> pd.DataFrame:
    df = df.copy()

    return df[(df["N"]) == 1].reset_index(drop=True)


def load_config(path: str) -> dict:
    """Load YAML configuration file."""
    with open(path, "r") as yaml_file:
        config = yaml.safe_load(yaml_file)
    return config

def get_num_sequences(fasta_file):
    """Count sequences in a FASTA file."""
    return sum(1 for _ in SeqIO.parse(fasta_file, "fasta"))



def run_hmmsearch(input_file, output_dir, hmm_file):
    """
    Run hmmsearch using a given HMM file against a FASTA database and save output.
    """
    hmm_file = Path(hmm_file)
    fasta_file = Path(input_file)
    hmmer_output_file = Path(os.path.join(output_dir, "hmmer2ODD_hits.out"))
    hmmer_output_csv = Path(os.path.join(output_dir, "hmmer2ODD_hits.csv"))


    cmd = ["hmmsearch", str(hmm_file), str(fasta_file)]

    print("Running command:", " ".join(cmd), ">", str(hmmer_output_file))
    
    # Open the output file and redirect stdout there
    with hmmer_output_file.open("w") as out_f:
        result = subprocess.run(cmd, stdout=out_f, stderr=subprocess.PIPE, text=True)
    
    if result.returncode != 0:
        print("hmmsearch failed!")
        print("STDERR:", result.stderr)
        raise subprocess.CalledProcessError(result.returncode, cmd, result.stderr)
    
    print(f"hmmsearch completed. Output saved to {hmmer_output_file}")

    df = parse_hmmsearch_output(hmmer_output_file)
    df.to_csv(hmmer_output_csv)
    if len(df) > 0:
        df = filter_hmmer_df(df)

    return df



def main(input_file, config):
    """Run the full filter_seq_sim pipeline for a single FASTA file."""
    output_dir = config["output_dir"]

    # if enabled, check if the output has already been calculated and skip the computation
    if config["reuse_existing"]:
        hmmer_hits = os.path.join(output_dir, "hmmer2ODD_hits.out")
        hmmer_hits_csv = os.path.join(output_dir, "hmmer2ODD_hits.csv")
        filtered_hmmer_fasta = os.path.join(output_dir, "filter_seq_sim_hmmer.fasta")
        if all(os.path.isfile(f) for f in [hmmer_hits, hmmer_hits_csv, filtered_hmmer_fasta]):
            print(f"SKIPPING HMMER COMPUTATION FOR {os.path.basename(output_dir)}: existing results found.")
            return 

    os.makedirs(output_dir, exist_ok=True)

    print(f"Processing {input_file}")
    print(f"  → input size: {get_num_sequences(input_file)} ")

    hmmer_cfg = config["tools"]["hmmer"]
    hmmer_output_df = run_hmmsearch(
        input_file=input_file,
        output_dir=output_dir,
        hmm_file=hmmer_cfg["domain_model"]
    )

    accs = get_accessions_safe(hmmer_output_df, "Sequence")
    filter_fasta_by_headers(input_file=input_file, output_dir=output_dir, accession_list=accs)
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



if __name__ == "__main__":
    import time

    config = load_config("/Users/michellealexander/projects/FALGs/configs/filter.yml")
    input_path = config["input"]

    start_time = time.time()
    print(f"Starting filter_seq_sim - hmmer process at {time.ctime(start_time)}")

    if os.path.isfile(input_path):
        main(input_path, config)
    elif os.path.isdir(input_path):
        process_directory(input_path, config)
    else:
        raise ValueError("Input must be a FASTA file or a directory containing FASTA files.")

    end_time = time.time()
    print(f"Finished in {end_time - start_time:.1f} seconds")

