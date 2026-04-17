from __future__ import annotations

from pathlib import Path
import subprocess

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment


def run_mafft(
    input_fasta: str | Path, output_fasta: str | Path, threads: int | None = None
) -> None:
    """Run MAFFT on ``input_fasta`` and write the alignment to ``output_fasta``.

    This is a thin wrapper around the ``mafft`` CLI and assumes it is
    installed and on the PATH, as in the original implementation.
    """

    input_fasta = Path(input_fasta)
    output_fasta = Path(output_fasta)

    cmd = ["mafft", "--auto"]
    if threads is not None:
        if threads < 1:
            raise ValueError(f"threads must be >= 1, got {threads}")
        cmd.extend(["--thread", str(int(threads))])
    cmd.append(str(input_fasta))
    with output_fasta.open("w") as out_f:
        result = subprocess.run(cmd, stdout=out_f, stderr=subprocess.PIPE, text=True)

    if result.returncode != 0:
        raise RuntimeError(f"mafft failed for {input_fasta}: {result.stderr}")


def trim_msa_by_gap_fraction(
    input_alignment: str | Path,
    output_alignment: str | Path,
    gap_threshold: float = 0.9,
) -> None:
    """Trim an alignment by removing columns with too many gaps.

    Any column whose gap fraction exceeds ``gap_threshold`` is removed.
    """

    input_alignment = Path(input_alignment)
    output_alignment = Path(output_alignment)

    alignment = AlignIO.read(str(input_alignment), "fasta")
    n_seqs = len(alignment)
    if n_seqs == 0:
        output_alignment.touch()
        return

    keep_columns: list[int] = []
    aln_len = alignment.get_alignment_length()

    for idx in range(aln_len):
        column = alignment[:, idx]
        gap_fraction = column.count("-") / n_seqs
        if gap_fraction <= gap_threshold:
            keep_columns.append(idx)

    if not keep_columns:
        # If all columns were filtered out, fall back to the original
        AlignIO.write(alignment, str(output_alignment), "fasta")
        return

    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    trimmed_records: list[SeqRecord] = []
    for row in alignment:
        seq = "".join(row.seq[i] for i in keep_columns)
        trimmed_records.append(SeqRecord(Seq(seq), id=row.id, description=""))

    trimmed_alignment = MultipleSeqAlignment(trimmed_records)
    AlignIO.write(trimmed_alignment, str(output_alignment), "fasta")
