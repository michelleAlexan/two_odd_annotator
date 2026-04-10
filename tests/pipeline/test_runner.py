#%%
import pytest
from pathlib import Path
from two_odd_annotator.pipeline.runner import Runner

# under tests/data/ there are two test files
test_data_path = Path(__file__).parents[1] / "data"
# save the test run restults under tests/results
test_results_path = Path(__file__).parents[1] / ".results"


pipeline = Runner(
    input_path=test_data_path,
    output_base_dir=test_results_path,
    reuse_existing=False,  # force re-run of all steps for testing purposes
)
pipeline.run()



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
        seq_sim_method="diamond",
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
        assert (sub / "diamond_results.tsv").exists()
        assert (sub / "filtered_diamond_hits.csv").exists()
        assert (sub / "filtered_diamond.fasta").exists()
