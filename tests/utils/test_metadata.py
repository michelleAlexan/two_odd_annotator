import pytest
import yaml
from two_odd_annotator.utils.metadata import infer_species_from_dir_name, create_metadata


def test_infer_species_valid():
    assert infer_species_from_dir_name("Solanum_tuberosum") == "Solanum tuberosum"
    assert infer_species_from_dir_name("Arabidopsis_thaliana") == "Arabidopsis thaliana"


def test_infer_species_with_path():
    assert (
        infer_species_from_dir_name("/some/path/Solanum_tuberosum")
        == "Solanum tuberosum"
    )


def test_infer_species_invalid():
    with pytest.raises(ValueError):
        infer_species_from_dir_name("Solanum")



def test_create_metadata_valid_tax_id(tmp_path):
    dir_name = "Arabidopsis_thaliana"
    dir_path = tmp_path / dir_name
    dir_path.mkdir()

    create_metadata(str(dir_path))

    meta_file = dir_path / "metadata.yml"
    assert meta_file.is_file()

    # Load YAML properly
    data = yaml.safe_load(meta_file.read_text())

    # Check species
    assert data["species"] == "Arabidopsis thaliana"

    # Check tax_id dict
    assert data["tax_id"] == {"Arabidopsis thaliana": 3702}

def test_create_metadata_invalid_tax_id(tmp_path):
    dir_name = "Unknown_species"
    dir_path = tmp_path / dir_name
    dir_path.mkdir()

    with pytest.raises(ValueError) as e:
        create_metadata(str(dir_path))


