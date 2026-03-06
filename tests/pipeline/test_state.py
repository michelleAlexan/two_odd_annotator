import pytest

from two_odd_annotator.pipeline.state import State
from two_odd_annotator.constants import (
    DIAMOND_RESULTS,
    FILTERED_DIAMOND_HITS,
    FILTERED_DIAMOND_FASTA,
    HMMER_RESULTS,
    FILTERED_HMMER_HITS,
    FILTERED_HMMER_FASTA,
    BLASTP_RESULTS,
    FILTERED_BLASTP_HITS,
    FILTERED_BLASTP_FASTA,
)


# ------------ TESTS FOR _validate_fasta_input_exists() ---------------
def test__validate_fasta_input_exists_no_output_dir_fasta_in_input_fasta_paths(tmp_path):
    # Create a dummy input FASTA file
    dummy_fasta = tmp_path / "Zea_mays.fasta"
    dummy_fasta.write_text(">seq1\nACGT\n")

    # Create an output directory that is a subdirectory of the input directory
    output_dir = tmp_path / "results"
    output_dir.mkdir()

    # Initialize State with the input path and output directory
    state = State(input_path=tmp_path, output_base_dir=output_dir)

    # Assert that the input FASTA file is correctly identified and not confused with any files in the output directory
    assert state.input_fasta_paths == [dummy_fasta]

def test__validate_fasta_input_exists_input_path_equal_output_base_dir(tmp_path):
    with pytest.raises(ValueError, match= (f"Input path and output directory must be different: " \
                f"{tmp_path}")):
        State(input_path=tmp_path, output_base_dir=tmp_path)

def test__validate_fasta_input_exists_single_file_valid(tmp_path):
    file = tmp_path / "Zea_mays.fasta"
    file.write_text(">seq1\nACGT\n")
    state = State(input_path=file, output_base_dir=tmp_path / "results")
    assert state.input_fasta_paths == [file]

def test__validate_fasta_input_exists_single_file_non_existent(tmp_path):
    fake_file = tmp_path / "does_not_exist.pep.fasta"
    with pytest.raises(FileNotFoundError):
        State(input_path=fake_file, output_base_dir=tmp_path / "results")

def test__validate_fasta_input_exists_single_file_invalid_extension(tmp_path):
    txt_file = tmp_path / "file.pep.txt"
    txt_file.write_text("Hello\n")
    with pytest.raises(ValueError, match="is not a FASTA file"):
        State(input_path=txt_file, output_base_dir=tmp_path / "results")

def test__validate_fasta_input_exists_dir_with_fasta(tmp_path):
    fasta_file = tmp_path / "Zea_mays.fasta"
    fasta_file.write_text(">seq1\nACGT\n")
    
    state = State(input_path=tmp_path, output_base_dir=tmp_path / "results")
    assert state.input_fasta_paths == [fasta_file]

def test__validate_fasta_input_exists_dir_empty(tmp_path):
    with pytest.raises(ValueError, match="No FASTA files found"):
        State(input_path=tmp_path, output_base_dir=tmp_path / "results")

def test__validate_fasta_input_exists_dir_not_exist(tmp_path):
    fake_path = tmp_path / "does_not_exist"
    with pytest.raises(FileNotFoundError):
        State(input_path=fake_path, output_base_dir=tmp_path / "results")

def test__validate_fasta_input_exists_dir_recursive(tmp_path):
    sub1 = tmp_path / "subdir1"
    sub2 = tmp_path / "subdir2"
    sub1.mkdir()
    sub2.mkdir()
    f1 = sub1 / "Zea_mays.fasta"
    f2 = sub2 / "Solanum_tuberosum.fasta"
    f1.write_text(">seq1\nACGT\n")
    f2.write_text(">seq2\nACGT\n")

    state = State(input_path=tmp_path, output_base_dir=tmp_path / "results")
    assert f1 in state.input_fasta_paths
    assert f2 in state.input_fasta_paths
    assert len(state.input_fasta_paths) == 2


def test__validate_input_dir_duplicate_species(tmp_path):
    f1 = tmp_path / "Solanum_tuberosum.fasta"
    f2 = tmp_path / "Solanum_tuberosum.fa"
    f1.write_text(">seq1\nACGT\n")
    f2.write_text(">seq2\nACGT\n")

    with pytest.raises(ValueError, match="Duplicate species names found"):
        State(input_path=tmp_path, output_base_dir=tmp_path / "results")



# ------------ TESTS FOR _infer_species_from_file_name() ---------------
def test__infer_species_from_file_name(tmp_path):
    dummy_file = tmp_path / "Zea_mays.fasta"
    dummy_file.write_text(">seq1\nACGT\n")
    state = State(input_path=dummy_file, output_base_dir=tmp_path)

    assert state._infer_species_from_file_name("Arabidopsis_thaliana") == "Arabidopsis thaliana"
    assert state._infer_species_from_file_name("Echinochloa_crus-galli") == "Echinochloa crus-galli"
    assert state._infer_species_from_file_name("Oryza_sativa") == "Oryza sativa"

    with pytest.raises(ValueError, match="Could not infer species name"):
        state._infer_species_from_file_name("invalidfilename.fasta")


# ------------ TESTS FOR _get_base_name()  ---------------

def test__get_base_name(tmp_path):
    dummy_file = tmp_path / "Zea_mays.fasta"
    dummy_file.write_text(">seq1\nACGT\n")
    state = State(input_path=dummy_file, output_base_dir=tmp_path)

    assert state._get_base_name("Arabidopsis_thaliana.fasta") == "Arabidopsis_thaliana"
    assert state._get_base_name("Echinochloa_crus-galli.fa") == "Echinochloa_crus-galli"
    assert state._get_base_name("Oryza_sativa.pep.faa") == "Oryza_sativa"
    assert state._get_base_name("complex.name.with.dots.fasta") == "complex"
    assert state._get_base_name("no_extension") == "no_extension"


#------------ TESTS FOR _fetch_completed_pipeline_steps_for_subdir() ---------------

def test__fetch_completed_pipeline_steps_for_subdir_no_steps_completed(tmp_path):
    dummy_input_fasta = tmp_path / "Zea_mays.fasta"
    dummy_input_fasta.write_text(">seq1\nACGT\n")
    dummy_output_subdir= tmp_path / "Zea_mays"
    dummy_output_subdir.mkdir()

    state = State(input_path=dummy_input_fasta, output_base_dir=tmp_path)

    assert state.results == {
        "Zea_mays": {
            "seq_sim_filter": None,
            "annotate": None,
            "visualize": None
        }
    }

def test__fetch_completed_pipeline_steps_for_subdir_filter_completed(tmp_path):
    dummy_input_fasta = tmp_path / "Zea_mays.fasta"
    dummy_input_fasta.write_text(">seq1\nACGT\n")
    dummy_output_subdir= tmp_path / "Zea_mays"
    dummy_output_subdir.mkdir()

    (dummy_output_subdir / DIAMOND_RESULTS).touch()
    (dummy_output_subdir / FILTERED_DIAMOND_HITS).touch()
    (dummy_output_subdir / FILTERED_DIAMOND_FASTA).touch()

    state = State(input_path=dummy_input_fasta, output_base_dir=tmp_path)
    assert state.results["Zea_mays"]["seq_sim_filter"] == ["diamond"]


def test__fetch_completed_pipeline_steps_for_subdir_only_single_subdir_completed(tmp_path):
    # seq_sim_filter step completed for Zea_mays, but not for Solanum_tuberosum
    dummy_input_fasta = tmp_path / "Zea_mays.fasta"
    dummy_input_fasta.write_text(">seq1\nACGT\n")
    dummy_input_fasta2 = tmp_path / "Solanum_tuberosum.fasta"
    dummy_input_fasta2.write_text(">seq2\nACGT\n")

    output_dir = tmp_path / "results"
    output_dir.mkdir()
    dummy_output_subdir= output_dir / "Zea_mays"
    dummy_output_subdir.mkdir()

    (dummy_output_subdir / DIAMOND_RESULTS).touch()
    (dummy_output_subdir / FILTERED_DIAMOND_HITS).touch()
    (dummy_output_subdir / FILTERED_DIAMOND_FASTA).touch()  


    state = State(input_path=tmp_path, output_base_dir=output_dir)
    assert state.results == {
        "Zea_mays": {
            "seq_sim_filter": ["diamond"],
            "annotate": None,
            "visualize": None
        },
        "Solanum_tuberosum": {
            "seq_sim_filter": None,
            "annotate": None,
            "visualize": None
        }
    }   


