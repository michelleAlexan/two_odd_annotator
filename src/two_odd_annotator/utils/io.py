from pathlib import Path
from typing import Any

import yaml

METADATA_FILENAME = "metadata.yml"


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
        output_path = output_path / METADATA_FILENAME
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
        meta_path = path / METADATA_FILENAME
    elif path.is_file() and path.name == METADATA_FILENAME:
        meta_path = path
    else:
        raise FileNotFoundError(
            f"Expected metadata file '{METADATA_FILENAME}' in directory or as a file."
        )

    if not meta_path.is_file():
        raise FileNotFoundError(f"Metadata file not found: {meta_path}")

    with meta_path.open("r") as f:
        try:
            data = yaml.safe_load(f)
        except yaml.YAMLError as e:
            raise ValueError(f"Invalid YAML in metadata file: {meta_path}") from e

    return data or {}