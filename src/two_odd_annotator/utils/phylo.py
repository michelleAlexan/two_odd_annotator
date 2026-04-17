from __future__ import annotations

from pathlib import Path
import subprocess
import os
import shutil


def run_fasttree(
    input_alignment: str | Path, output_tree: str | Path, threads: int | None = None
) -> None:
    """Run FastTree on an alignment file to produce a Newick tree.

    This is a thin wrapper around the ``fasttree`` CLI and assumes it is
    installed and on the PATH, mirroring the behaviour of the original
    helper that lived in ``bio_tools``.
    """

    input_alignment = Path(input_alignment)
    output_tree = Path(output_tree)

    fasttree_bin = (
        shutil.which("FastTreeMP")
        or shutil.which("FastTree")
        or shutil.which("fasttree")
        or "fasttree"
    )
    cmd = [fasttree_bin, str(input_alignment)]
    env = None
    if threads is not None:
        if threads < 1:
            raise ValueError(f"threads must be >= 1, got {threads}")
        env = os.environ.copy()
        env["OMP_NUM_THREADS"] = str(int(threads))
    with output_tree.open("w") as out_f:
        result = subprocess.run(
            cmd, stdout=out_f, stderr=subprocess.PIPE, text=True, env=env
        )

    if result.returncode != 0:
        raise RuntimeError(f"FastTree failed for {input_alignment}: {result.stderr}")
