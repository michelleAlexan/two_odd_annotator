#%%
import argparse
import os
import subprocess
import timeit
from pathlib import Path
from sys import argv
from typing import Literal

import pandas as pd
from Bio import SeqIO

from bio_tools.files.fasta import clean_fasta_file
from two_odd_annotator.utils.io import load_config
from two_odd_annotator.utils.logging import init_log, log_line
from two_odd_annotator.utils.metadata import init_subdir

DEFAULT_CONFIG_PATH = Path(__file__).parents[2] / "configs" / "filter.yml"

BLAST_COLUMNS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore"
]


#%%

# =============================================================================
# Shared utilities
# =============================================================================

BLAST_COLUMNS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore"
]

def count_fasta_sequences(fasta_file: str) -> int:
    return sum(1 for _ in SeqIO.parse(fasta_file, "fasta"))


def write_filtered_fasta(input_fasta: str, output_fasta: str, ids: set):
    with open(output_fasta, "w") as out:
        for record in SeqIO.parse(input_fasta, "fasta"):
            if record.id in ids:
                SeqIO.write(record, out, "fasta")


# =============================================================================
#  DIAMOND OR BLAST 
# =============================================================================
def run_alignment(
    tool: Literal["diamond", "blastp"],
    input_file: str,
    output_dir: str,
    config: dict,
) -> pd.DataFrame:
    """
    Run either DIAMOND or BLASTP for a single input FASTA file
    and return the filtered results as a DataFrame.
    """

    os.makedirs(output_dir, exist_ok=True)

    # -----------------------
    # Load configuration
    # -----------------------
    reference_db = config["filter_tools"][tool]["reference_db"]

    threads = config["parameters"]["threads"]
    thresholds = config["parameters"]["thresholds_alignment"]

    output_tsv = os.path.join(output_dir, f"{tool}_results.tsv")
    filtered_csv = os.path.join(output_dir, f"filtered_{tool}_hits.csv")
    filtered_fasta = os.path.join(output_dir, f"filtered_{tool}.fasta")

    # -----------------------
    # Reuse existing results
    # -----------------------
    if config.get("reuse_existing") and all(
        os.path.isfile(f) for f in [output_tsv, filtered_csv, filtered_fasta]
    ):
        print(f"SKIPPING {tool.upper()} for {os.path.basename(output_dir)}")
        return pd.read_csv(filtered_csv)

    print(f"Processing {input_file}")
    print(f"  → input size: {count_fasta_sequences(input_file)}")

    # -----------------------
    # Build command
    # -----------------------
    if tool == "blastp":
        cmd = (
            f"blastp "
            f"-query {input_file} "
            f"-db {reference_db} "
            f"-out {output_tsv} "
            f"-evalue {thresholds['evalue']} "
            f"-num_threads {threads} "
            f"-outfmt '6 {' '.join(BLAST_COLUMNS)}'"
        )

    elif tool == "diamond":
        cmd = (
            f"diamond blastp "
            f"--query {input_file} "
            f"--db {reference_db} "
            f"--out {output_tsv} "
            f"--evalue {thresholds['evalue']} "
            f"--threads {threads} "
            f"--sensitive "
            f"--outfmt 6 {' '.join(BLAST_COLUMNS)}"
        )

    print("Run:", cmd)
    result = subprocess.run(cmd, shell=True)
    if result.returncode != 0:
        raise RuntimeError(f"{tool} failed")

    # -----------------------
    # Load results
    # -----------------------
    df = pd.read_csv(output_tsv, sep="\t", header=None, names=BLAST_COLUMNS)

    # -----------------------
    # Apply filtering
    # -----------------------
    df_filtered = (
        df[(df["pident"] >= thresholds["pident"]) &
           (df["length"] >= thresholds["length"]) &
           (df["bitscore"] >= thresholds["bitscore"])]
        .sort_values("bitscore", ascending=False)
        .groupby("qseqid")
        .head(thresholds["num_hits"])
        .reset_index(drop=True)
    )

    df_filtered.to_csv(filtered_csv, index=False)

    passed_ids = set(df_filtered["qseqid"])    # Save filtered CSV (even if empty)
    df_filtered.to_csv(filtered_csv, index=False)

    # Write filtered FASTA (empty if no hits)
    passed_ids = set(df_filtered["qseqid"])
    if passed_ids:
        write_filtered_fasta(input_file, filtered_fasta, passed_ids)
    else:
        # create empty FASTA if no sequences passed thresholds
        open(filtered_fasta, "w").close()


    print(f"  → {len(passed_ids)} sequences retained")

    return df_filtered



# =============================================================================
# HMMER 
# =============================================================================


def parse_hmmsearch_output(hmmout_path: str) -> pd.DataFrame:
    """
    Parse HMMER 'Scores for complete sequences' table into a DataFrame.
    Assumes table has at least 9 columns:
    full_Evalue, full_score, full_bias, bestdom_Evalue, bestdom_score, bestdom_bias, exp, N, Sequence, [Description]
    """
    rows = []
    in_table = False
    with open(hmmout_path) as f:
        for line in f:
            if "Scores for complete sequences" in line:
                in_table = False
            if line.strip().startswith("E-value"):
                in_table = True
                continue
            if "---" in line and in_table:
                continue
            if in_table and line.strip() == "":
                break
            if in_table:
                parts = line.strip().split()
                if len(parts) < 9:
                    continue
                rows.append({
                    "full_Evalue": float(parts[0]),
                    "full_score": float(parts[1]),
                    "full_bias": float(parts[2]),
                    "bestdom_Evalue": float(parts[3]),
                    "bestdom_score": float(parts[4]),
                    "bestdom_bias": float(parts[5]),
                    "exp": float(parts[6]),
                    "N": int(parts[7]),
                    "Sequence": parts[8],
                    "Description": " ".join(parts[9:]) if len(parts) > 9 else "",
                })

    df= pd.DataFrame(rows)
    
    return df


def run_hmmer(input_file: str, output_dir: str, config: dict) -> pd.DataFrame:
    """
    Run hmmsearch, filter using thresholds from config, and write filtered FASTA and CSV.
    """
    hmmer_cfg = config["filter_tools"]["hmmer"]
    thresholds = {key: float(val) for key, val in config["parameters"]["thresholds_hmmer"].items()}

    os.makedirs(output_dir, exist_ok=True)
    output_txt = Path(output_dir) / "hmmer_results.out"
    filtered_csv = Path(output_dir) / "filtered_hmmer_hits.csv"
    filtered_fasta = Path(output_dir) / "filtered_hmmer.fasta"

    # Skip if reuse_existing
    if config.get("reuse_existing") and all(f.exists() for f in [output_txt, filtered_csv, filtered_fasta]):
        print(f"Skipping HMMER: existing results found for {input_file}")
        return pd.read_csv(filtered_csv)

    # Run hmmsearch
    cmd = ["hmmsearch", hmmer_cfg["domain_model"], input_file]
    print("Run:", " ".join(cmd))
    with open(output_txt, "w") as out_f:
        result = subprocess.run(cmd, stdout=out_f, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"hmmsearch failed for {input_file}:\n{result.stderr}")

    # Parse results
    df = parse_hmmsearch_output(output_txt)
    if df.empty:
        # If no hits found, create empty CSV and FASTA
        pd.DataFrame(columns=[
            "full_Evalue","full_score","full_bias",
            "bestdom_Evalue","bestdom_score","bestdom_bias",
            "exp","N","Sequence","Description"
        ]).to_csv(filtered_csv, index=False)

        open(filtered_fasta, "w").close()  # create empty fasta
        print("  → No sequences passed HMMER filtering. Empty files created.")
        return df

    # Apply thresholds
    df_filtered = df[
        (df["full_Evalue"] <= thresholds["full_Evalue"]) &
        (df["bestdom_Evalue"] <= thresholds["bestdom_Evalue"]) &
        (df["full_score"] >= thresholds["full_score"]) &
        (df["bestdom_score"] >= thresholds["bestdom_score"]) &
        (df["N"] == thresholds["N"])
    ].reset_index(drop=True)

    # Save CSV
    df_filtered.to_csv(filtered_csv, index=False)

    # Filter FASTA
    keep_ids = set(df_filtered["Sequence"])
    if not df.empty:
        write_filtered_fasta(input_file, filtered_fasta, keep_ids)

    print(f"HMMER: {len(keep_ids)} sequences passed thresholds")

    return df_filtered



# =============================================================================
# Pipeline logic 
# =============================================================================

def run_single_file(input_path: Path|str, output_base_dir: Path|str, config: dict, method: str, logfile_path: str) -> None:
    """ Run the full filtering pipeline for a single input FASTA file.
    This is the core function that processes one FASTA file, runs DIAMOND or HMMER or BLASTP filtering, and handles metadata updates. 
    It is called by the main CLI function for each file found in the input path.
    """

    input_path = Path(input_path)
    output_base_dir = Path(output_base_dir)

    species, tax_id, subdir = init_subdir(input_path, output_base_dir)

    log_line(logfile_path, f"Processing {input_path.name}")

    clean_fasta_file(input_path, species=species, tax_info=tax_id, output_dir=subdir)

    if method in ("blastp", "diamond"):
        run_alignment(method, str(input_path), str(subdir), config)

    elif method == "hmmer":
        run_hmmer(str(input_path), str(subdir), config)

    log_line(logfile_path, f"Finished {input_path.name}")


def run(input_path, output_dir, config_path, method):

    input_path = Path(input_path)
    output_dir = Path(output_dir)
    config = load_config(config_path)

    logfile = init_log(str(output_dir), "filter.log")

    start = timeit.default_timer()

    if input_path.is_file():
        run_single_file(input_path, output_dir, config, method, logfile)

    else:
        fasta_files = [
            os.path.join(root, f)
            for root, _, files in os.walk(input_path)
            for f in files
            if f.lower().endswith((".fasta", ".fa", ".faa")) and ".cds" not in f
        ]

        for f in fasta_files:
            run_single_file(f, output_dir, config, method, logfile)

    elapsed = (timeit.default_timer() - start) / 60
    log_line(logfile, f"Total time: {elapsed:.2f} minutes")


# =============================================================================
# CLI
# =============================================================================

def main_cli(argv=None):
    """CLI entry point for the filter pipeline.

    This orchestrates DIAMOND or BLASTP or HMMER filtering for one or many
    FASTA files and handles logging and tax ID mapping.
    """

    parser = argparse.ArgumentParser(
        description="Filter plant protein FASTA files for 2ODD candidates."
    )

    parser.add_argument("--input-path", required=True, help="Path to an input FASTA file or a directory containing plant protein FASTA files.")
    parser.add_argument("--output-dir", required=True, help="Base output directory for all results.")
    parser.add_argument("--config-path", required=True, default=DEFAULT_CONFIG_PATH, help="Path to YAML configuration file (thresholds, DBs, etc.)")
    parser.add_argument(
        "--method",
        choices=["diamond", "blastp", "hmmer"],
        default="diamond",
        help="Which filtering method to run.",
    )

    args = parser.parse_args(argv)

    run(
        input_path=args.input_path,
        output_dir=args.output_dir,
        config_path=args.config_path,
        method=args.method,
    )


if __name__ == "__main__":
    # if the script is run from the command line, execute the main CLI function with the provided arguments.
    if len(argv) > 1:
        main_cli(argv[1:])
    else:
        # If no arguments are provided, run with default settings for testing or development.
        print("Running with default settings.")
        run(
            input_path="/Users/michellealexander/projects/two_odd_annotator/tests/data",
            output_dir="/Users/michellealexander/projects/two_odd_annotator/tests/results",
            config_path=DEFAULT_CONFIG_PATH,
            method="hmmer"
        )

 


