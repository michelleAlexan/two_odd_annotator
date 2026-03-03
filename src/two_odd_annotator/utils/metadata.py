
#%%
from pathlib import Path
from datetime import datetime, timezone
from typing import Optional, Dict, Any
from bio_tools.taxa.taxonomy import map_scientific_notation_to_tax_id

from two_odd_annotator.utils.io import write_metadata


METADATA_FILENAME = "metadata.yml"


def _now_iso() -> str:
    """Return current UTC time as ISO 8601 string."""

    return datetime.now(timezone.utc).isoformat()


def infer_species_from_file_name(file_name: Path) -> str:
    """Infer species name from an input FASTA file name.

    Assumes the file is named according to the Latin species name, e.g.
    "Solanum_tuberosum.fasta" or
    "Solanum_tuberosum.pep.fasta"
    -> "Solanum tuberosum".
    """

    file_name = Path(file_name)

    stem = file_name.name
    while True:
        stem_path = Path(stem)
        if stem_path.suffix == "":
            break
        stem = stem_path.stem

    # Replace underscores with spaces
    name = stem.replace("_", " ").strip()

    parts = name.split()
    if len(parts) < 2:
        raise ValueError(f"Could not infer species name from file: {file_name}")

    return name


def init_subdir(input_fasta_file: Path, output_base_dir: Path) -> None:
    """it is assumed that the input fasta file is named according to the Latin species name.
    For downstream processing, the tax id needs to be inferred from the filename and stored in the metadata."""
    inferred_species = infer_species_from_file_name(input_fasta_file)
    tax_id_dict = map_scientific_notation_to_tax_id(inferred_species, raise_on_error=True)
    inferred_species, tax_id = list(tax_id_dict.keys())[0], list(tax_id_dict.values())[0]
    metadata = {
        "creation_timestamp": _now_iso(),
        "species": inferred_species,
        "tax_id": tax_id,
    }
    # create a subdirectory named as the scientific name of the species
    subdir = output_base_dir / inferred_species.replace(" ", "_")
    subdir.mkdir(parents=True, exist_ok=True)

    write_metadata(output_path=subdir / METADATA_FILENAME, metadata=metadata)
    return inferred_species, tax_id, subdir
    















def ensure_base_metadata(output_dir: str, input_file: str) -> Dict[str, Any]:
    """Ensure a metadata file exists with creation/species/tax_id.

    If metadata.yml already exists in the output directory, it is loaded
    and returned unchanged. Otherwise a new metadata structure is
    created, written, and returned.
    """

    metadata = load_metadata(output_dir)
    if metadata:
        return metadata

    species = infer_species_from_filename(input_file)
    tax_id = resolve_tax_id(species)

    metadata = {
        "creation_timestamp": _now_iso(),
        "species": species,
        "tax_id": tax_id,
    }
    save_metadata(output_dir, metadata)
    return metadata


def update_metadata(output_dir: str, updates: Dict[str, Any]) -> Dict[str, Any]:
    """Merge updates into existing metadata and save.

    Shallow-merge only; nested dicts will be overwritten.
    Returns the updated metadata dict.
    """

    metadata = load_metadata(output_dir)
    metadata.update(updates)
    save_metadata(output_dir, metadata)
    return metadata
