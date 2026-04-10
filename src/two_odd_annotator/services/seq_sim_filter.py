# %%
from logging import config
import os
import subprocess
from pathlib import Path
from typing import Literal

import pandas as pd
from Bio import SeqIO

BLAST_COLUMNS = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
]

# %%
# =============================================================================
# Helper functions
# =============================================================================


def count_fasta_sequences(fasta_file: str) -> int:
    return sum(1 for _ in SeqIO.parse(fasta_file, "fasta"))


def write_filtered_fasta(
    input_fasta: str, output_fasta: str, ids: set, seq_len_thresh: int | None
) -> None:
    median_2ODD_seq_len = 349
    with open(output_fasta, "w") as out:
        for record in SeqIO.parse(input_fasta, "fasta"):
            if record.id in ids:
                if seq_len_thresh is not None:
                    if (
                        len(record.seq) < median_2ODD_seq_len - seq_len_thresh
                        or len(record.seq) > median_2ODD_seq_len + seq_len_thresh
                    ):
                        continue
                SeqIO.write(record, out, "fasta")


# =============================================================================
#  ALIGNMENT FUNCTION (DIAMOND OR BLAST)
# =============================================================================


def run_alignment(
    tool: Literal["diamond", "blastp"],
    input_file: str,
    output_dir: str,
    config: dict,
) -> pd.DataFrame:
    """
    Run either DIAMOND or BLASTP for a single input FASTA file
    and return the pre-filtered results as a DataFrame.
    """

    os.makedirs(output_dir, exist_ok=True)

    # Load configuration file
    reference_db = config["filter_tools"][tool]["reference_db"]

    threads = config["parameters"]["threads"]
    thresholds = config["parameters"]["thresholds_alignment"]
    seq_len_thresh = config["pipeline"]["seq_len_thresh"]
    if seq_len_thresh == -1:
        seq_len_thresh = None

    output_tsv = os.path.join(output_dir, f"{tool}_results.tsv")
    filtered_csv = os.path.join(output_dir, f"filtered_{tool}_hits.csv")
    filtered_fasta = os.path.join(output_dir, f"filtered_{tool}.fasta")

    # Reuse existing results
    if config.get("reuse_existing") and all(
        os.path.isfile(f) for f in [output_tsv, filtered_csv, filtered_fasta]
    ):
        print(f"SKIPPING {tool.upper()} for {os.path.basename(output_dir)}")
        return pd.read_csv(filtered_csv)

    print(f"Processing {input_file}")
    print(f"  → input size: {count_fasta_sequences(input_file)}")

    # Build command
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

    # Load results
    df = pd.read_csv(output_tsv, sep="\t", header=None, names=BLAST_COLUMNS)

    # Apply filtering
    df_filtered = (
        df[
            (df["pident"] >= thresholds["pident"])
            & (df["length"] >= thresholds["length"])
            & (df["bitscore"] >= thresholds["bitscore"])
        ]
        .sort_values("bitscore", ascending=False)
        .groupby("qseqid")
        .head(thresholds["num_hits"])
        .reset_index(drop=True)
    )

    passed_ids = set(df_filtered["qseqid"])
    df_filtered.to_csv(filtered_csv, index=False)

    # Write filtered FASTA (empty if no hits)
    passed_ids = set(df_filtered["qseqid"])
    if passed_ids:
        write_filtered_fasta(
            input_file, filtered_fasta, passed_ids, seq_len_thresh=seq_len_thresh
        )
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
                rows.append(
                    {
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
                    }
                )

    df = pd.DataFrame(rows)

    return df


def run_hmmer(input_file: str, output_dir: str, config: dict) -> pd.DataFrame:
    """
    Run hmmsearch, filter using thresholds from config, and write filtered FASTA and CSV.
    """
    hmmer_cfg = config["filter_tools"]["hmmer"]
    thresholds = {
        key: float(val) for key, val in config["parameters"]["thresholds_hmmer"].items()
    }
    seq_len_thresh = config["pipeline"]["seq_len_thresh"]
    if seq_len_thresh == -1:
        seq_len_thresh = None

    os.makedirs(output_dir, exist_ok=True)
    output_txt = Path(output_dir) / "hmmer_results.out"
    filtered_csv = Path(output_dir) / "filtered_hmmer_hits.csv"
    filtered_fasta = Path(output_dir) / "filtered_hmmer.fasta"

    # Skip if reuse_existing
    if config.get("reuse_existing") and all(
        f.exists() for f in [output_txt, filtered_csv, filtered_fasta]
    ):
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
        pd.DataFrame(
            columns=[
                "full_Evalue",
                "full_score",
                "full_bias",
                "bestdom_Evalue",
                "bestdom_score",
                "bestdom_bias",
                "exp",
                "N",
                "Sequence",
                "Description",
            ]
        ).to_csv(filtered_csv, index=False)

        open(filtered_fasta, "w").close()  # create empty fasta
        print("  → No sequences passed HMMER filtering. Empty files created.")
        return df

    # Apply thresholds
    df_filtered = df[
        (df["full_Evalue"] <= thresholds["full_Evalue"])
        & (df["bestdom_Evalue"] <= thresholds["bestdom_Evalue"])
        & (df["full_score"] >= thresholds["full_score"])
        & (df["bestdom_score"] >= thresholds["bestdom_score"])
        & (df["N"] >= thresholds["N"])
    ].reset_index(drop=True)

    # Save CSV
    df_filtered.to_csv(filtered_csv, index=False)

    # Filter FASTA
    keep_ids = set(df_filtered["Sequence"])
    if not df.empty:
        write_filtered_fasta(
            input_file, filtered_fasta, keep_ids, seq_len_thresh=seq_len_thresh
        )

    print(f"HMMER: {len(keep_ids)} sequences passed thresholds")

    return df_filtered


# ============================================================================
# SEQUENCE SIMILIRATY FILTERING (Pipeline logic)
# =============================================================================
def run(
    input_path: Path | str, subdir: Path | str, config: dict, seq_sim_method: str
) -> None:
    """Run the full sequence similarity filtering pipeline for a single input FASTA file.
    This is the core function that processes one FASTA file, runs DIAMOND or HMMER or BLASTP filtering, and handles metadata updates.
    It is called by the main CLI function for each file found in the input path.
    """

    input_path = Path(input_path)
    subdir = Path(subdir)

    if seq_sim_method in ("blastp", "diamond"):
        run_alignment(seq_sim_method, str(input_path), str(subdir), config)

    elif seq_sim_method == "hmmer":
        run_hmmer(str(input_path), str(subdir), config)
