from pathlib import Path
from typing import Any

import json

from Bio import SeqIO

from two_odd_annotator.constants import METADATA_YML, CLEAN_FASTA_HEADERS_JSON

import yaml


def load_config(path: str | Path) -> dict[str, Any]:
    """Load a YAML configuration file."""
    path = Path(path)

    if not path.is_file():
        raise FileNotFoundError(f"Config file not found: {path}")

    with path.open("r") as f:
        try:
            data = yaml.safe_load(f)
        except yaml.YAMLError as e:
            raise ValueError(f"Invalid YAML in config file: {path}") from e

    return data or {}


def write_metadata(output_path: Path, metadata: dict[str, Any]) -> Path:
    """Write metadata dict as YAML into the given path.

    Returns
    -------
    Path
        Path to the written metadata file.
    """
    if output_path.is_dir():
        output_path = output_path / METADATA_YML
    else:
        output_path = Path(output_path)

    with output_path.open("w") as f:
        yaml.safe_dump(metadata, f, sort_keys=False)

    return output_path


def load_metadata(path: str | Path) -> dict[str, Any]:
    """Load metadata YAML from a directory or file.

    Raises
    ------
    FileNotFoundError
        If metadata file is missing.
    ValueError
        If YAML is invalid.
    """
    path = Path(path)

    if path.is_dir():
        meta_path = path / METADATA_YML
    elif path.is_file() and path.name == METADATA_YML:
        meta_path = path
    else:
        raise FileNotFoundError(
            f"Expected metadata file '{METADATA_YML}' in directory or as a file."
        )

    if not meta_path.is_file():
        raise FileNotFoundError(f"Metadata file not found: {meta_path}")

    with meta_path.open("r") as f:
        try:
            data = yaml.safe_load(f)
        except yaml.YAMLError as e:
            raise ValueError(f"Invalid YAML in metadata file: {meta_path}") from e

    return data or {}


def _build_clean_id(orig_header: str, scientific_sp_name: str, tax_id: int | str) -> str:
    """Construct a cleaned FASTA identifier embedding species and taxid.

    Heuristics:
    - If the header already ends with ``__<taxid>``, keep it as-is.
    - If it matches ``<acc> <gene> [Species name]``, produce
      ``<acc>_<gene>_Species_name_with_underscores__taxid``.
    - Otherwise, take the first whitespace-separated token and append
      ``__taxid``.
    """

    header = orig_header.strip()
    tax_str = str(tax_id)

    if header.endswith(f"__{tax_str}"):
        return header

    import re

    m = re.match(r"^(?P<acc>\S+)\s+(?P<gene>\S+)\s+\[(?P<species>[^\]]+)\]", header)
    if m:
        acc = m.group("acc")
        gene = m.group("gene")
        species = m.group("species").replace(" ", "_")
        return f"{acc}_{gene}_{species}__{tax_str}"

    base = header.split()[0]
    return f"{base}__{tax_str}"


def write_clean_fasta_with_taxid(
    input_fasta_path: str | Path,
    output_dir: str | Path,
    output_fasta_name: str,
    scientific_sp_name: str,
    tax_info: int | str,
) -> Path:
    """Copy a FASTA file while normalising headers and recording a mapping.

    For each record in ``input_fasta_path`` this function:
    - Derives a new identifier embedding the species name and tax ID.
    - Writes the updated records into ``<output_dir>/<output_fasta_name>.fasta``.
    - Writes a JSON mapping of original headers to cleaned IDs and tax info
      into ``clean_fasta_headers.json`` in ``output_dir``.
    """

    input_fasta_path = Path(input_fasta_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    tax_id = int(tax_info)

    output_fasta_path = output_dir / f"{output_fasta_name}.fasta"
    mapping_path = output_dir / CLEAN_FASTA_HEADERS_JSON

    mapping: dict[str, dict[str, Any]] = {}
    records = []

    for record in SeqIO.parse(str(input_fasta_path), "fasta"):
        orig_header = record.description
        clean_id = _build_clean_id(orig_header, scientific_sp_name, tax_id)

        record.id = clean_id
        record.description = ""
        records.append(record)

        mapping[orig_header] = {
            "clean_id": clean_id,
            "tax_id": tax_id,
            "scientific_name": scientific_sp_name,
        }

    if records:
        SeqIO.write(records, str(output_fasta_path), "fasta")
    else:
        output_fasta_path.touch()

    with mapping_path.open("w") as f:
        json.dump(mapping, f, indent=2)

    return output_fasta_path
