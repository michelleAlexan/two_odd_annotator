import pytest
from pathlib import Path

from two_odd_annotator.services.annotate import run
from two_odd_annotator.utils.io import load_config
from two_odd_annotator.constants import (
    DEFAULT_CONFIG_PATH,
    ANNOTATION_FASTA,
    ANNOTATION_MSA,
    ANNOTATION_MSA_TRIM,
    ANNOTATION_TREE
)

RESULTS_DIR = Path(__file__).parents[1] /  "results" 

# delete existing annotation results files before running tests
if (RESULTS_DIR / ANNOTATION_FASTA).exists():
    (RESULTS_DIR / ANNOTATION_FASTA).unlink()
if (RESULTS_DIR / ANNOTATION_MSA).exists():
    (RESULTS_DIR / ANNOTATION_MSA).unlink()
if (RESULTS_DIR / ANNOTATION_MSA_TRIM).exists():
    (RESULTS_DIR / ANNOTATION_MSA_TRIM).unlink()
if (RESULTS_DIR / ANNOTATION_TREE).exists():
    (RESULTS_DIR / ANNOTATION_TREE).unlink()

config = load_config(Path(__file__).parents[2] / DEFAULT_CONFIG_PATH)
config["annotate"]["ingroup"] = Path(__file__).parents[2] / "data" / "2ODDs" / "characterized_2ODDs.fasta"

def test_annotate():
    run(
        result_dir=RESULTS_DIR, 
        config=config, 
        seq_sim_method="hmmer"
    )

    # check that annotation results files were created
    assert (RESULTS_DIR / ANNOTATION_FASTA).exists(), "Annotation FASTA file was not created."
    assert (RESULTS_DIR / ANNOTATION_MSA).exists(), "Annotation MSA file was not created."
    assert (RESULTS_DIR / ANNOTATION_MSA_TRIM).exists(), "Annotation trimmed MSA file was not created."
    assert (RESULTS_DIR / ANNOTATION_TREE).exists(), "Annotation tree file was not created."


