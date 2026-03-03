import pytest
import yaml
from two_odd_annotator.utils.metadata import infer_species_from_file_name, init_subdir


def test_infer_species_valid():
    assert infer_species_from_file_name("Solanum_tuberosum.pep.fasta") == "Solanum tuberosum"
    assert infer_species_from_file_name("Arabidopsis_thaliana.fasta") == "Arabidopsis thaliana"


def test_infer_species_with_path():
    assert (
        infer_species_from_file_name("/some/path/Solanum_tuberosum.fasta")
        == "Solanum tuberosum"
    )


def test_infer_species_invalid():
    with pytest.raises(ValueError):
        infer_species_from_file_name("Solanum.fasta")


def test_init_subdir_valid_tax_id(tmp_path):
    input_path = "some/input/folder/Arabidopsis_thaliana.pep.fasta"
    output_base_path = tmp_path 
    init_subdir(input_path, output_base_path)
    subdir = output_base_path / "Arabidopsis_thaliana"

    meta_file = subdir / "metadata.yml"
    assert meta_file.is_file()

    # Load YAML properly
    data = yaml.safe_load(meta_file.read_text())

    # Check species
    assert data["species"] == "Arabidopsis thaliana"

    # Check tax_id dict
    assert data["tax_id"] == 3702

def test_init_subdir_invalid_tax_id(tmp_path):
    dir_name = "Unknown_species"
    dir_path = tmp_path / dir_name
    dir_path.mkdir()

    with pytest.raises(ValueError) as e:
        init_subdir(str(dir_path), tmp_path)


