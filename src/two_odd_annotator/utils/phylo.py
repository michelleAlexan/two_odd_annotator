from __future__ import annotations

from pathlib import Path
import subprocess


def run_fasttree(input_alignment: str | Path, output_tree: str | Path) -> None:
    """Run FastTree on an alignment file to produce a Newick tree.

    This is a thin wrapper around the ``fasttree`` CLI and assumes it is
    installed and on the PATH, mirroring the behaviour of the original
    helper that lived in ``bio_tools``.
    """

    input_alignment = Path(input_alignment)
    output_tree = Path(output_tree)

    cmd = ["fasttree", str(input_alignment)]
    with output_tree.open("w") as out_f:
        result = subprocess.run(cmd, stdout=out_f, stderr=subprocess.PIPE, text=True)

    if result.returncode != 0:
        raise RuntimeError(f"FastTree failed for {input_alignment}: {result.stderr}")
