
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


def infer_species_from_dir_name(dir_name: str) -> str:
    """Infer species name from an input directory name.

    Assumes the directory is named according to the Latin species name, e.g.
    "Solanum_tuberosum" -> "Solanum tuberosum".
    """

    stem = Path(dir_name).stem
    # Replace underscores with spaces and title-case words
    name = stem.replace("_", " ").strip()

    if len(name.split()) < 2:
        raise ValueError(f"Could not infer species name from directory name: {dir_name}")

    return name


def create_metadata(dir_path: str) -> None:
    """Write the scientific species name and the corresponding tax_id as a metadata YAML into the given directory."""
    
    inferred_species = infer_species_from_dir_name(dir_path)
    tax_id = map_scientific_notation_to_tax_id(inferred_species, raise_on_error=True)
    metadata = {
        "creation_timestamp": _now_iso(),
        "species": inferred_species,
        "tax_id": tax_id,
    }
    write_metadata(dir_path, metadata)
    















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
