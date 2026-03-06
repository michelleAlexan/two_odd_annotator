import pytest
from pathlib import Path

from two_odd_annotator.utils.io import write_metadata, load_metadata


def test_write_and_load_metadata_roundtrip(tmp_path: Path):
    metadata = {
        "Arabidopsis_thaliana": {
            "scientific_sp_name": "Arabidopsis thaliana",
            "tax_id": 3702,
            "original_file_path": "/path/to/arabidopsis.fasta"
        },
    }


    write_metadata(tmp_path, metadata)
    loaded = load_metadata(tmp_path)

    assert loaded == metadata

def test_load_metadata_missing(tmp_path: Path):
    with pytest.raises(FileNotFoundError):
        load_metadata(tmp_path)
    

