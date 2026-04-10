from pathlib import Path
import pytest

from two_odd_annotator.pipeline.runner import Runner
from two_odd_annotator.utils.io import load_config
from two_odd_annotator.constants import DEFAULT_CONFIG_PATH

TEST_INPUT_DIR = Path(__file__).parents[1] / "data" 

test_config_path = Path(__file__).parents[1] / "config" / "test_config.yml"


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
>NP_001190266.1_FLS_Arabidopsis_thaliana__3702
MEVERVQDISSSSLLTEAIPLEFIRSEKEQPAITTFRGPTPAIPVVDLSDPDEESVRRAV
VKASEEWGLFQVVNHGIPTELIRRLQDVGRKFFELPSSEKESVAKPEDSKDIEGYGTKLQ
KDPEGKKAWVDHLFHRIWPPSCVNYRFWPKNPPEYREVNEEYAVHVKKLSETLLGILSDG
LGLKRDALKEGLGGEMAEYMMKINYYPPCPRPDLALGVPAHTDLSGITLLVPNEVPGLQV
FKDDHWFDAEYIPSAVIVHIGDQILRLSNGRYKNVLHRTTVDKEKTRMSWPVFLEPPREK
IVGPLPELTGDDNPPKFKPFAFKDYSYRKLNKLPLD
"""

expected_echinochloa_filtered_2ODD_fasta = ""

def test_all_filter_methods(tmp_path):
    output_dir = tmp_path / "results"
    input_dir = TEST_INPUT_DIR 



    for method in ["diamond", "blastp", "hmmer"]:
        pipline = Runner(
            input_path=input_dir,
            output_base_dir=output_dir,
            seq_sim_method=method,
            reuse_existing=False,  # force re-run of all steps for testing purposes
            config_path=test_config_path,
            step="filter_seq_sim"

        )

        pipline.run()

        # check that all filtered FASTA files contain the expected headers and sequences for each method
        with open(output_dir / "Arabidopsis_thaliana" / f"filtered_{method}.fasta") as f:
            filtered_fasta = f.read()
            assert filtered_fasta == expected_arabidposis_filtered_2ODD_fasta, f"Filtered FASTA does not match expected for Arabidopsis thaliana with method {method}"
        with open(output_dir / "Echinochloa_crus-galli" / f"filtered_{method}.fasta") as f:
            filtered_fasta = f.read()
            assert filtered_fasta == expected_echinochloa_filtered_2ODD_fasta, f"Filtered FASTA does not match expected for Echinochloa crus-galli with method {method}"



def test_seq_len_thresh_default_100(tmp_path):
    # fasta content with seq length of 186
    test_fasta_path = Path(__file__).parents[1] / "data" / "Vitis_rotundifolia.pep.fasta"
    pipline = Runner(
        input_path=test_fasta_path,
        output_base_dir=tmp_path / "results", 
        seq_sim_method="diamond",
        config_path=test_config_path,
        step="filter_seq_sim"
    )
    pipline.run()

    # check that all filtered FASTA files contain the expected headers and sequences for each method
    with open(tmp_path / "results" / "Vitis_rotundifolia" / f"filtered_diamond.fasta") as f:
        filtered_fasta = f.read()
        assert filtered_fasta == "", "Filtered FASTA should be empty since sequence length is outside of default threshold of 100 from median 2ODD sequence length of 349"

EXPECTED_FASTA = """>lcl_CM058159.1_cds_KAJ9692488.1_17292__103349
MAVYDKDRLGVHNTVVGNQNIIDKPYHGYTLERSVASLHQGLGIDNCKKLIQSFAKVIFP
PLIPTKSVMSILHQNQVNGLEIETKDGKWIGYEPLTPSLFSQGSKIHPPKHQVIMKGNEA
RYSLRLFSFTKGLIKIPKELVDDHHPLQFQSFDHIDLLISFAQKKVEILRVLLKPTEALL
KISLI*
"""

def test_seq_len_thresh_override(tmp_path):
    # fasta content with seq length of 186
    test_fasta_path = Path(__file__).parents[1] / "data" / "Vitis_rotundifolia.pep.fasta"
    pipline = Runner(
        input_path=test_fasta_path,
        output_base_dir=tmp_path / "results",
        seq_sim_method="diamond",
        seq_len_thresh=200,
        config_path=test_config_path,
        step="filter_seq_sim"

    )
    pipline.run()

    # check that all filtered FASTA files contain the expected headers and sequences for each method
    with open(tmp_path / "results" / "Vitis_rotundifolia" / f"filtered_diamond.fasta") as f:
        filtered_fasta = f.read()
        assert filtered_fasta == EXPECTED_FASTA, "Filtered FASTA should contain the sequence since it is within the overridden threshold of 200 from median 2ODD sequence length of 349"


def test_seq_len_thresh_no_filtering(tmp_path):
    # fasta content with seq length of 186
    test_fasta_path = Path(__file__).parents[1] / "data" / "Vitis_rotundifolia.pep.fasta"
    pipline = Runner(
        input_path=test_fasta_path,
        output_base_dir=tmp_path / "results",
        seq_sim_method="diamond",
        seq_len_thresh=-1,
        step="filter_seq_sim",
        config_path=test_config_path
    )
    pipline.run()

    # check that all filtered FASTA files contain the expected headers and sequences for each method
    with open(tmp_path / "results" / "Vitis_rotundifolia" / f"filtered_diamond.fasta") as f:
        filtered_fasta = f.read()
        assert filtered_fasta == EXPECTED_FASTA, "Filtered FASTA should contain the sequence since no length filtering should be applied when threshold is set to -1"

