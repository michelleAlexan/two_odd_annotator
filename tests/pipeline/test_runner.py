#%%
import pytest
from pathlib import Path
from two_odd_annotator.pipeline.runner import Runner

from two_odd_annotator.constants import (
    DIAMOND_RESULTS ,
    HMMER_RESULTS,
    BLASTP_RESULTS,

    FILTERED_DIAMOND_HITS ,
    FILTERED_HMMER_HITS,
    FILTERED_BLASTP_HITS, 

    FILTERED_DIAMOND_FASTA,
    FILTERED_HMMER_FASTA,
    FILTERED_BLASTP_FASTA,


    ANNOTATION_CSV,
    ANNOTATION_FASTA,
    ANNOTATION_MSA,
    ANNOTATION_MSA_TRIM,
    ANNOTATION_TREE,
)

# under tests/data/ there are two test files
test_data_path = Path(__file__).parents[1] / "data"
# save the test run restults under tests/results
test_results_path = Path(__file__).parents[1] / ".results"

test_config_path = Path(__file__).parents[1] / "config" / "test_config.yml"



# ---------------- TEST EXPECTED FOLDER STRUCTURE ----------------

# The expected folder structure after running the filter pipeline with diamond method 
# on the test datasets is as follows:
#
# results/
# ├── Arabidopsis_thaliana/
# │   ├── metadata.yml
# │   ├── clean_Arabidopsis_thaliana.pep.fasta
# │   ├── clean_fasta_headers.json
# │   ├── diamond_results.tsv
# │   ├── filtered_diamond_hits.csv
# │   ├── filtered_diamond.fasta
# ├── Echinochloa_crus-galli/
# │   ├── metadata.yml
# │   ├── clean_Echinochloa_crus-galli.pep.fasta
# │   ├── clean_fasta_headers.json
# │   ├── diamond_results.tsv
# │   ├── filtered_diamond_hits.csv
# │   ├── filtered_diamond.fasta


# %%

def test_runner_pipeline(tmp_path):

    output_dir = tmp_path / "results"

    runner = Runner(
        input_path=test_data_path,
        output_base_dir=output_dir,
        seq_sim_method="all",  # run all methods for testing purposes
        config_path=test_config_path
    )

    runner.run()

    # ---- Top-level checks ----
    assert output_dir.exists()
    # assert (output_dir / "log.log").exists()

    # ---- Subdirectory checks ----
    expected_species = [
        "Arabidopsis_thaliana",
        "Echinochloa_crus-galli",
    ]

    for sp in expected_species:
        sub = output_dir / sp

        assert sub.exists()
        assert (sub / "metadata.yml").exists()
        assert (sub / f"clean_{sp}.fasta").exists()
        assert (sub / "clean_fasta_headers.json").exists()
        assert (sub / DIAMOND_RESULTS).exists()
        assert (sub / FILTERED_DIAMOND_HITS).exists()
        assert (sub / FILTERED_DIAMOND_FASTA).exists()
        assert (sub / HMMER_RESULTS).exists()
        assert (sub / FILTERED_HMMER_HITS).exists()
        assert (sub / FILTERED_HMMER_FASTA).exists()
        assert (sub / BLASTP_RESULTS).exists()
        assert (sub / FILTERED_BLASTP_HITS).exists()
        assert (sub / FILTERED_BLASTP_FASTA).exists()

    # expected annotation results files
    assert (output_dir / ANNOTATION_CSV).exists()
    assert (output_dir / ANNOTATION_FASTA).exists()
    assert (output_dir / ANNOTATION_MSA).exists()
    assert (output_dir / ANNOTATION_MSA_TRIM).exists()
    assert (output_dir / ANNOTATION_TREE).exists()
