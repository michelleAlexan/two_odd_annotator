from pathlib import Path
import pytest

from two_odd_annotator.pipeline.runner import Runner
from two_odd_annotator.utils.io import load_config
from two_odd_annotator.constants import DEFAULT_CONFIG_PATH

TEST_INPUT_DIR = Path(__file__).parents[1] / "data" 



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
        )

        pipline.run()

        # check that all filtered FASTA files contain the expected headers and sequences for each method
        with open(output_dir / "Arabidopsis_thaliana" / f"filtered_{method}.fasta") as f:
            filtered_fasta = f.read()
            assert filtered_fasta == expected_arabidposis_filtered_2ODD_fasta, f"Filtered FASTA does not match expected for Arabidopsis thaliana with method {method}"
        with open(output_dir / "Echinochloa_crus-galli" / f"filtered_{method}.fasta") as f:
            filtered_fasta = f.read()
            assert filtered_fasta == expected_echinochloa_filtered_2ODD_fasta, f"Filtered FASTA does not match expected for Echinochloa crus-galli with method {method}"

