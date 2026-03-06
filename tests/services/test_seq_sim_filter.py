from pathlib import Path
import pytest

from two_odd_annotator.services.seq_sim_filter import run

TEST_INPUT_DIR = Path(__file__).parent / "data" 
TEST_OUTPUT_DIR = Path(__file__).parent / "results"
TEST_CONFIG = Path(__file__).parent / "data" / "filter.yml"
print(f"Using test config: {TEST_CONFIG}")  


expected_arabidposis_filtered_2ODD_fasta = """>sp|Q96323.1_ANS_Arabidopsis_thaliana__3702
MVAVERVESLAKSGIISIPKEYIRPKEELESINDVFLEEKKEDGPQVPTIDLKNIESDDE
KIRENCIEELKKASLDWGVMHLINHGIPADLMERVKKAGEEFFSLSVEEKEKYANDQATG
KIQGYGSKLANNASGQLEWEDYFFHLAYPEEKRDLSIWPKTPSDYIEATSEYAKCLRLLA
TKVFKALSVGLGLEPDRLEKEVGGLEELLLQMKINYYPKCPQPELALGVEAHTDVSALTF
ILHNMVPGLQLFYEGKWVTAKCVPDSIVMHIGDTLEILSNGKYKSILHRGLVNKEKVRIS
WAVFCEPPKDKIVLKPLPEMVSVESPAKFPPRTFAQHIEHKLFGKEQEELVSEKND
>NP_190692.1_F3H_Arabidopsis_thaliana__3702
MAPGTLTELAGESKLNSKFVRDEDERPKVAYNVFSDEIPVISLAGIDDVDGKRGEICRQI
VEACENWGIFQVVDHGVDTNLVADMTRLARDFFALPPEDKLRFDMSGGKKGGFIVSSHLQ
GEAVQDWREIVTYFSYPVRNRDYSRWPDKPEGWVKVTEEYSERLMSLACKLLEVLSEAMG
LEKESLTNACVDMDQKIVVNYYPKCPQPDLTLGLKRHTDPGTITLLLQDQVGGLQATRDN
GKTWITVQPVEGAFVVNLGDHGHFLSNGRFKNADHQAVVNSNSSRLSIATFQNPAPDATV
YPLKVREGEKAILEEPITFAEMYKRKMGRDLELARLKKLAKEERDHKEVDKPVDQIFA
>NP_199094.1_DFR_Arabidopsis_thaliana__3702
>NP_001190266.1_FLS_Arabidopsis_thaliana__3702
MEVERVQDISSSSLLTEAIPLEFIRSEKEQPAITTFRGPTPAIPVVDLSDPDEESVRRAV
VKASEEWGLFQVVNHGIPTELIRRLQDVGRKFFELPSSEKESVAKPEDSKDIEGYGTKLQ
KDPEGKKAWVDHLFHRIWPPSCVNYRFWPKNPPEYREVNEEYAVHVKKLSETLLGILSDG
LGLKRDALKEGLGGEMAEYMMKINYYPPCPRPDLALGVPAHTDLSGITLLVPNEVPGLQV
FKDDHWFDAEYIPSAVIVHIGDQILRLSNGRYKNVLHRTTVDKEKTRMSWPVFLEPPREK
IVGPLPELTGDDNPPKFKPFAFKDYSYRKLNKLPLD
"""

expected_echinochloa_filtered_2ODD_fasta = ""

######################## TESTS FOR SUCCESSFUL RUNS OF THE FILTER PIPELINE ########################
#
#
# ---------------- TEST EXPECTED FOLDER STRUCTURE ----------------

# The expected folder structure after running the filter pipeline with diamond method 
# on the test datasets is as follows:
#
# results/
# ├── filter.log
# ├── Arabidopsis_thaliana/
# │   ├── metadata.yaml
# │   ├── diamond_results.tsv
# │   ├── filtered_diamond_hits.csv
# │   ├── filtered_diamond.fasta
# ├── Echinochloa_crus-galli/
# │   ├── metadata.yaml
# │   ├── diamond_results.tsv
# │   ├── filtered_diamond_hits.csv
# │   ├── filtered_diamond.fasta



def test_filter_pipeline_diamond(tmp_path):
    output_dir = tmp_path / "results"

    run(
        str(TEST_INPUT_DIR),
        str(output_dir),
        str(TEST_CONFIG),
        method="diamond",
    )

    # ---- Top-level checks ----
    assert output_dir.exists()
    assert (output_dir / "filter.log").exists()

    # ---- Subdirectory checks ----
    expected_species = [
        "Arabidopsis_thaliana",
        "Echinochloa_crus-galli",
    ]

    for sp in expected_species:
        sub = output_dir / sp

        assert sub.exists()
        assert (sub / "metadata.yml").exists()
        assert (sub / "diamond_results.tsv").exists()
        assert (sub / "filtered_diamond_hits.csv").exists()
        assert (sub / "filtered_diamond.fasta").exists()
    

def test_filter_pipeline_hmmer(tmp_path):
    output_dir = tmp_path / "results"

    run(
        str(TEST_INPUT_DIR),
        str(output_dir),
        str(TEST_CONFIG),
        method="hmmer",
    )

    # ---- Top-level checks ----
    assert output_dir.exists()
    assert (output_dir / "filter.log").exists()

    # ---- Subdirectory checks ----
    expected_species = [
        "Arabidopsis_thaliana",
        "Echinochloa_crus-galli",
    ]

    for sp in expected_species:
        sub = output_dir / sp

        assert sub.exists()
        assert (sub / "metadata.yml").exists()
        assert (sub / "hmmer_results.out").exists()
        assert (sub / "filtered_hmmer_hits.csv").exists()
        assert (sub / "filtered_hmmer.fasta").exists()


def test_filter_pipeline_blastp(tmp_path):
    output_dir = tmp_path / "results"

    run(
        str(TEST_INPUT_DIR),
        str(output_dir),
        str(TEST_CONFIG),
        method="blastp",
    )

    # ---- Top-level checks ----
    assert output_dir.exists()
    assert (output_dir / "filter.log").exists()

    # ---- Subdirectory checks ----
    expected_species = [
        "Arabidopsis_thaliana",
        "Echinochloa_crus-galli",
    ]

    for sp in expected_species:
        sub = output_dir / sp

        assert sub.exists()
        assert (sub / "metadata.yml").exists()
        assert (sub / "blastp_results.tsv").exists()
        assert (sub / "filtered_blastp_hits.csv").exists()
        assert (sub / "filtered_blastp.fasta").exists()


def test_all_methods_pipeline(tmp_path):
    """
    Run the full filter pipeline for all three methods sequentially
    and check that the results folder contains the correct structure.
    """

    output_dir = tmp_path / "results"

    # Run each method in sequence
    for method in ["diamond", "blastp", "hmmer"]:
        run(
            str(TEST_INPUT_DIR),
            str(output_dir),
            str(TEST_CONFIG),
            method=method
        )

    # ---- Top-level checks ----
    assert output_dir.exists(), "Results folder was not created"
    assert (output_dir / "filter.log").exists(), "Log file missing"

    # ---- Species subdirectories ----
    expected_species = ["Arabidopsis_thaliana", "Echinochloa_crus-galli"]
    for sp in expected_species:
        subdir = output_dir / sp
        assert subdir.exists(), f"Subdirectory {sp} missing"
        assert (subdir / "metadata.yml").exists(), f"Metadata missing in {sp}"

        # For each method, check that files exist (or empty but present)
        for method in ["diamond", "blastp", "hmmer"]:
            if method == "diamond" or method == "blastp":
                expected_result_file = f"{method}_results.tsv"
                expected_csv = f"filtered_{method}_hits.csv"
                expected_fasta = f"filtered_{method}.fasta"
            elif method == "hmmer":
                expected_result_file = f"{method}_results.out"
                expected_csv = f"filtered_{method}_hits.csv"
                expected_fasta = f"filtered_{method}.fasta"
            


            assert (subdir / expected_result_file).exists(), f"{method} results file missing in {sp}"
            assert (subdir / expected_csv).exists(), f"{method} filtered CSV missing in {sp}"
            assert (subdir / expected_fasta).exists(), f"{method} filtered FASTA missing in {sp}"




"""
######################## TESTS FOR UNSUCCESSFUL RUNS OF THE FILTER PIPELINE ########################

def test_run_filter():
    # Test the run function with a small test dataset"
    input_dir = mock_input_files["input_dir"]
    output_dir = RESULTS_DIR
    output_dir.mkdir(parents=True, exist_ok=True)

    # Run the filter pipeline
    run(str(input_dir), str(output_dir), str(TEST_CONFIG), method=method, plots=False)

    # check that loggning file is created
    log_file = output_dir / "filter.log"
    assert log_file.exists(), "Log file not created"

    # Check that output files are created for each input FASTA file
    expected_subdirs = ["Arabidopsis_thaliana", "Echinochloa_crus-galli"]
    for subdir in expected_subdirs:
        subdir_path = output_dir / subdir
        assert subdir_path.exists() and subdir_path.is_dir(), f"Subdirectory {subdir} not created"

        # Check for presence of expected output files (e.g. BLASTP results)
        blast_result_file = subdir_path / "blast_results.tsv"
        filtered_hits = subdir_path / "filtered_blast_hits.csv"
        filtered_fasta = subdir_path / "filtered_blast.fasta"
        assert blast_result_file.exists(), f"BLAST results not found in {subdir}"
        assert filtered_hits.exists(), f"Filtered BLAST hits not found in {subdir}"
        assert filtered_fasta.exists(), f"Filtered BLAST FASTA not found in {subdir}"

        # Check for presence of expected output files (e.g. HMMER results)
        hmmer_hits = subdir_path / "hmmer_results.out"
        hmmer_hits_csv = subdir_path / "filtered_hmmer_hits.csv"
        filtered_hmmer_fasta = subdir_path / "filtered_hmmer.fasta"
        assert hmmer_hits.exists(), f"HMMER hits not found in {subdir}"
        assert hmmer_hits_csv.exists(), f"HMMER hits CSV not found in {subdir_path}"
        assert filtered_hmmer_fasta.exists(), f"Filtered HMMER FASTA not found in {subdir_path}"


    expected_headers_in_filtered_blast_fasta_arabidposis = {
        ">spQ96323.1_ANS_Arabidopsis_thaliana__3702", 
        ">NP_190692.1_F3H_Arabidopsis_thaliana__3702", 
        ">NP_001190266.1_FLS_Arabidopsis_thaliana__3702"
    }

    expected_headers_in_filtered_blast_fasta_echinochloa = {}

    # Check that the filtered BLAST FASTA files contain the expected headers
    with open(output_dir / "Arabidopsis_thaliana" / "filtered_blast.fasta") as f:
        headers = {line.strip() for line in f if line.startswith(">")}
        assert headers == expected_headers_in_filtered_blast_fasta_arabidposis, "Filtered BLAST FASTA headers do not match expected for Arabidopsis thaliana"
    with open(output_dir / "Echinochloa_crus-galli" / "filtered_blast.fasta") as f:
        headers = {line.strip() for line in f if line.startswith(">")}
        assert headers == expected_headers_in_filtered_blast_fasta_echinochloa, "Filtered BLAST FASTA headers do not match expected for Echinochloa crus-galli"
"""