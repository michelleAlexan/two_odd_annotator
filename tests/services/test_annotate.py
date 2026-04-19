#%%
import json

import pytest

import pandas as pd
from pathlib import Path
from ete4 import PhyloTree


from two_odd_annotator.utils.io import load_config
from two_odd_annotator.constants import (
    DEFAULT_CONFIG_PATH,
    ANNOTATION_FASTA,
    ANNOTATION_MSA,
    ANNOTATION_MSA_TRIM,
    ANNOTATION_TREE, 
)

from two_odd_annotator.services.annotate import (
    reverse_major_minor_2ODD_dict,
    split_seqs_by_2ODD_membership,
    assign_2ODD_props,
    assign_plant_group_props,
    create_annotation_fasta, 
    from_fasta_to_nwk,
    build_distance_lookup, 
    get_clusters, 
    compute_cluster_neighbors,
    two_odd_id_to_cluster_indices, 
    seq_to_cluster_idx,
    build_parent_map,
    cluster_meta_info,
    get_candidate_to_char_baits_df
)


def test_run_reuse_existing_derives_candidate_headers(tmp_path, monkeypatch):
    """Regression test for reruns with reuse_existing=True.

    When `completed_annotation_steps["annotation_fasta"]` is True, `run()` skips
    `create_annotation_fasta()`. It should still derive candidate/ingroup IDs from
    the existing annotation FASTA so downstream steps can proceed.
    """

    from two_odd_annotator.services import annotate as annotate_module

    class DummyPhyloTree:
        def __init__(self, _fh, sp_naming_function=None):
            self.sp_naming_function = sp_naming_function

        def annotate_ncbi_taxa(self, taxid_attr="species"):
            return {}, {}, {}

        def ladderize(self):
            return None

        def write(self, outfile=None, **kwargs):
            return None

    # Minimal ingroup + annotation FASTA
    ingroup_fasta = tmp_path / "ingroup.fasta"
    ingroup_fasta.write_text(">ING1\nMA\n")

    (tmp_path / ANNOTATION_FASTA).write_text(">ING1\nMA\n>CAND1\nMA\n")
    (tmp_path / ANNOTATION_TREE).write_text("(ING1:0.1,CAND1:0.1);\n")

    test_config = {
        "annotate": {
            "bait_sequence_collection": ingroup_fasta,
            "major_minor_2ODD_ids": major_minor_2ODDs_path,
            "major_2ODDs_functional_characterization": Path(__file__).parents[1]
            / "config"
            / "major_2ODD_char_info.json",
        },
        "pipeline": {"reuse_existing": True},
        "parameters": {"threads": 1},
    }

    completed_steps = {
        "annotation_fasta": True,
        "annotation_msa": True,
        "annotation_msa_trim": True,
        "annotation_tree": True,
    }

    captured: dict[str, set[str]] = {}

    def fake_assign_2ODD_props(tree, seq_to_2ODD_id, candidate_headers):
        captured["candidate_headers"] = set(candidate_headers)
        raise RuntimeError("STOP_AFTER_HEADERS")

    monkeypatch.setattr(annotate_module, "PhyloTree", DummyPhyloTree)
    monkeypatch.setattr(annotate_module, "assign_plant_group_props", lambda tree: None)
    monkeypatch.setattr(annotate_module, "assign_2ODD_props", fake_assign_2ODD_props)

    with pytest.raises(RuntimeError, match="STOP_AFTER_HEADERS"):
        annotate_module.run(
            result_dir=tmp_path,
            config=test_config,
            completed_annotation_steps=completed_steps,
        )

    assert captured["candidate_headers"] == {"CAND1"}


RESULTS_DIR = Path(__file__).parents[1] /  ".results" 

config = load_config(Path(__file__).parents[2] / DEFAULT_CONFIG_PATH)
config["annotate"]["bait_sequence_collection"] = Path(__file__).parents[2] / "data" / "2ODDs" / "characterized_2ODDs.fasta"

major_minor_2ODDs_path = Path(__file__).parents[1] /"config" / "major_minor_2ODD_ids_manual.json"
major_minor_2ODDs_dict = json.load(open(major_minor_2ODDs_path))
parent_map = build_parent_map(list(major_minor_2ODDs_dict["major_2ODDs"].keys()))

seq_to_2ODD_id = reverse_major_minor_2ODD_dict(major_minor_2ODD_dict=major_minor_2ODDs_dict)

def compare_nested_lists(result, expected):
    """
    Compare two lists of lists and print detailed differences.
    """
    print("Comparing nested lists...\n")
    max_len = max(len(result), len(expected))
    
    for i in range(max_len):
        try:
            r = result[i]
        except IndexError:
            print(f"Result missing list at index {i}: expected {expected[i]}")
            continue
        try:
            e = expected[i]
        except IndexError:
            print(f"Expected missing list at index {i}: result has {r}")
            continue
        
        if r != e:
            print(f"Difference at sublist index {i}:")
            print(f"  Result:   {r}")
            print(f"  Expected: {e}")
            # Check element-wise
            min_len = min(len(r), len(e))
            for j in range(min_len):
                if r[j] != e[j]:
                    print(f"    Element {j} differs: result={r[j]}, expected={e[j]}")
            # Extra elements
            if len(r) > len(e):
                print(f"    Extra elements in result: {r[min_len:]}")
            if len(e) > len(r):
                print(f"    Missing elements in result: {e[min_len:]}")
    print("\nComparison complete.")



#%%
def test_reverse_major_minor_2ODD_dict():
    input_dict = {
        "major_2ODDs": {
            "2ODD01": ["seq1", "seq2"],
            "2ODD02": ["seq3"]
        },
        "minor_2ODDs": {
            "2ODD_minor01": ["seq4", "seq5"]
        }
    }

    expected = {
        "seq1": "2ODD01",
        "seq2": "2ODD01",
        "seq3": "2ODD02",
        "seq4": "2ODD_minor01",
        "seq5": "2ODD_minor01"
    }

    result = reverse_major_minor_2ODD_dict(input_dict)

    assert result == expected


#%% Test file creations


def test_create_annotation_fasta():
    if (RESULTS_DIR / ANNOTATION_FASTA).exists():
        (RESULTS_DIR / ANNOTATION_FASTA).unlink()   

    
    create_annotation_fasta(
        results_dir=RESULTS_DIR,
        ingroup_2ODD_fasta=config["annotate"]["bait_sequence_collection"],
        output_fasta=RESULTS_DIR / ANNOTATION_FASTA,
        seq_sim_method="hmmer"
    )

    # check that annotation results files were created
    assert (RESULTS_DIR / ANNOTATION_FASTA).exists(), "Annotation FASTA file was not created."


def test_from_fasta_to_nwk():
    if (RESULTS_DIR / ANNOTATION_MSA).exists():
        (RESULTS_DIR / ANNOTATION_MSA).unlink()
    if (RESULTS_DIR / ANNOTATION_MSA_TRIM).exists():
        (RESULTS_DIR / ANNOTATION_MSA_TRIM).unlink()
    if (RESULTS_DIR / ANNOTATION_TREE).exists():
        (RESULTS_DIR / ANNOTATION_TREE).unlink()
    from_fasta_to_nwk(fasta_path=RESULTS_DIR / ANNOTATION_FASTA,
                      msa_path=RESULTS_DIR / ANNOTATION_MSA,
                        msa_trim_path=RESULTS_DIR / ANNOTATION_MSA_TRIM,
                        tree_path=RESULTS_DIR / ANNOTATION_TREE)
    

    assert (RESULTS_DIR / ANNOTATION_MSA).exists(), "Annotation MSA file was not created."
    assert (RESULTS_DIR / ANNOTATION_MSA_TRIM).exists(), "Annotation trimmed MSA file was not created."
    assert (RESULTS_DIR / ANNOTATION_TREE).exists(), "Annotation tree file was not created."


# %% test assign props functions

create_annotation_fasta(
    results_dir=RESULTS_DIR,
    ingroup_2ODD_fasta=config["annotate"]["bait_sequence_collection"],
    output_fasta=RESULTS_DIR / ANNOTATION_FASTA,
    seq_sim_method="hmmer"
    )

from_fasta_to_nwk(fasta_path=RESULTS_DIR / ANNOTATION_FASTA,
                    msa_path=RESULTS_DIR / ANNOTATION_MSA,
                    msa_trim_path=RESULTS_DIR / ANNOTATION_MSA_TRIM,
                    tree_path=RESULTS_DIR / ANNOTATION_TREE)

# load test tree
test_with_char_tree = RESULTS_DIR / ANNOTATION_TREE
t_char_2ODDs = PhyloTree(open(test_with_char_tree), sp_naming_function=lambda name: name.split('__')[-1])
tax2names, tax2lineages, tax2rank = t_char_2ODDs.annotate_ncbi_taxa(taxid_attr='species')

def test_assign_2ODD_props():
    candidate_headers, _ = split_seqs_by_2ODD_membership(seq_to_2ODD_id=seq_to_2ODD_id, tree=t_char_2ODDs)

    assign_2ODD_props(
        tree=t_char_2ODDs,
        seq_to_2ODD_id=seq_to_2ODD_id,
        candidate_headers=candidate_headers
    )

    assert t_char_2ODDs["PV023584__F3H__flavonoid_pathway__981085"].props["two_odd_id"] == "2ODD15"
    assert t_char_2ODDs["sp|Q96323.1_ANS_Arabidopsis_thaliana__3702"].props["two_odd_id"] == "candidate"

def test_assign_plant_group_props():
    assign_plant_group_props(tree=t_char_2ODDs)

    assert t_char_2ODDs["At4g10500__S3H__salicylic_acid_metabolism__3702"].props["plant_group"] == "Dicots"
    assert t_char_2ODDs["Os04g49210__S5H__salicylic_acid_metabolism__4530"].props["plant_group"] == "Monocots"


#%% ====================TEST CLUSTERING WITHOUT NESTED 2ODD IDs====================

test_anno_tree_path = Path(__file__).parents[1] / "data" / "test_anno_tree.nwk"
t = PhyloTree(open(test_anno_tree_path), sp_naming_function=lambda name: name.split('__')[-1])
tax2names, tax2lineages, tax2rank = t.annotate_ncbi_taxa(taxid_attr='species')
candidate_headers, ingroup_headers = split_seqs_by_2ODD_membership(seq_to_2ODD_id=seq_to_2ODD_id, tree=t)
assign_2ODD_props(
    tree=t, seq_to_2ODD_id=seq_to_2ODD_id, 
    candidate_headers=candidate_headers
)
assign_plant_group_props(tree=t)
# explore_2ODD_IDs(t)



expected_clusters_without_candidates = [
    [
        t["lcl_NC_084852.1_cds_XP_061960601.1_3168__3691"] # 2ODD38
    ], 
    [
        t["rna-XM_031626300.2__210225"], # 2ODD34
    ], 
    [
        t["CKAN_00595100__337451"], 
        t["Aristolochia_fimbriatalcl_CM034078.1_cds_KAG9459387.1_3895__158543"], 
        t["Gobar.D13G239800.1__3634"],
        t["Lj8A611G43.1__105884"], 
        t["FvH4_2g29990.t1__57918"],
        t['Casgl23S04257__3522'],
        t['CiLak.07G220400.1__32201'],
        t['lcl_NC_065570.1_cds_XP_050209022.1_2785__3986'],
        t['Lus10013130_PACid-23164132__4006'],
    ], 
    [
        t['evm_27.model.AmTr_v1.0_scaffold00174.4__13333'],
    ], 
    [
        t['RZC71586__3469'],
    ], 
    [
        t['lcl_NW_010729080.1_cds_XP_010249917.1_5454__4432'],
        t['rna-XM_042628350.1__60698'],
        t['lcl_CM056811.1_cds_KAJ8637692.1_11945__3435'],
        t['Aristolochia_fimbriatalcl_CM034078.1_cds_KAG9459386.1_3894__158543'],
    ], 
    [
        t['Kaladp1006s0012.1.v1.1__63787'],
        t['VIT_204s0008g04920.2__29760'],
        t['Ptrif.0001s0318.1__37690'],
        t['Atru_chr7_2342__47965'],
        t['lcl_NW_019168159.1_cds_XP_022728803.1_47147__66656'],
        t['lcl_NC_045134.1_cds_XP_031407529.1_33286__22663'],
        t['GWHTACBH000860__413952'],
        t['Potri.006G248000.1__3694'],
        t['lcl_NC_065570.1_cds_XP_050209042.1_2787__3986'],
        t['Lus10008097_PACid-23169790__4006'],
        t['Acc06446.1__3625'],
        t['Acc16585.1__3625'],
        t['pveT_jg24143.t1__170927'],
        t['pveT_jg28701.t1__170927'],
        t['lcl_NC_083383.1_cds_XP_048332159.2_17459__326968'],
        t['Cnepa31456.1__79760'],
        t['Umino40998.1__262084'],
        t['FvH4_2g29991.t1__57918'],
        t['Lj1g0022518.1__34305'],
        t['Lj4g0022664.1__34305'],
        t['Medtr1g011600.1__3880'],
        t['Ca_10933__3827'],
        t['Ler.1DRT.1g009880.1__41257'],
        t['Ler.1DRT.3g073930.1__41257'],
        t['Medtr3g108520.1__3880'],
        t['Ca_23071__3827'],
        t['Ca_01791__3827'],
        t['Medtr4g021360.1__3880'],
        t['Ca_05192__3827'],
        t['Medtr4g021380.1__3880'],
        t['lcl_OX451737.1_cds_CAI8600332.1_14923__3906']
    ], 
    [
        t["pveT_jg20218.t1__170927"], # minor 2ODD cluster
        t["pveT_jg20220.t1__170927"],
    ], 
    [
        t["Acc09929.1__3625"], # 2ODD34
    ], 
    [
        t["pveT_jg33784.t1__170927"] # minor 2ODD cluster
    ], 
    [
        t['Eucgr.C03724.1__71139'], # 2ODD38
        t['lcl_NC_045130.1_cds_XP_031390300.1_19564__22663'],
        t['Eucgr.C04147.1__71139'],
    ], 
    [
        t['Acc23206.1__3625'],
        t['lcl_CM027379.1_cds_KAF9612985.1_11655__261450'],
        t['RZC71605__3469'],
        t['rna-XM_042647017.1__60698'],
    ], 
    [
        t['lcl_CM027379.1_cds_KAF9612671.1_9820__261450']
    ],
    [
        t['Lj7A541T77.1__105884'] # minor 2ODD cluster
    ],
    [
        t['Pg_S7309.1__4054'],  # 2ODD34
        t['DCAR_015859__4039'],
    ], 
    [
        t['CsatW809539.1__3659'], 
    ], 
    [
        t['lcl_NW_026775580.1_cds_XP_059633894.1_39517__4283'],
        t['Dglom27812.1__34297'],
        t['lcl_CM059866.1_cds_KAK0595785.1_1317__4024'],
        t['Ptrif.0007s1329.1__37690'],
        t['Potri.001G176500.1__3694'],
        t['lcl_NC_065574.1_cds_XP_050234707.1_20560__3986'],
        t['Umino20162.1__262084'],
        t['lcl_NC_083379.1_cds_XP_015868728.2_4750__326968'],
        t['lcl_NC_084800.1_cds_XP_062086870.1_35676__3486'],
        t['lcl_NC_045131.1_cds_XP_031396687.1_21997__22663'],
        t['Eucgr.J03022.1__71139'],
        t['VIT_209s0002g05290.1__29760'],
        t['lcl_NW_026775580.1_cds_XP_059636056.1_39514__4283'],
        t['Acc08999.1__3625'],
        t['Acc07679.1__3625'],
        t['FSB010835101__28930'],
        t['DCAR_015858__4039']
    ], 
    [
        t['Casgl236S24615__3522'], # minor 2ODD cluster
        t['CiLak.15G093300.1__32201'],
        t['Qurub.10G093800.1__3512'],
        t['lcl_CM025849.1_cds_KAB1226463.1_1752__262757'],
        t['Tcord05348.2__703396'],
        t['HU01G00826.1__176265'],
        t['Sc13g0006470.01__3999'],
        t['TRINITY_DN108901_c0_g1_i1__223224'],
        t['CM008280.1.CM008280.1.g987.t1__62330']
    ], 
    [
        t['pveT_jg9723.t1__170927'] # 2ODD34
    ],
    [
        t["Cc06_g15940.1__49390"], # minor 2ODD cluster
    ], 
    [
        t["Gobar.D05G137100.1__3634"],
        t["lcl_OX459118.1_cds_CAI9089985.1_2909__43536"],
    ],
    [
        t['TRINITY_DN197356_c3_g1_i1__4102'], # 2ODD34
        t['lcl_CM061488.1_cds_KAK1430368.1_13046__13708'],
        t['Oeu000043.1__4146'],
        t['lcl_NW_026137546.1_cds_XP_051123007.1_5823__175694'],
        t['lcl_NW_025952766.1_cds_XP_047954044.1_46391__49212'],
        t['lcl_NC_028637.1_cds_XP_015064474.1_2338__28526'],
        t['Cc00_g25210.1__49390'],
        t['rna-gnl_WGS-JAMLDZ_EVM0036354.1__4058'],
        t['lcl_OX459121.1_cds_CAI9101309.1_14233__43536']
    ],
    [
        t["lcl_NC_031989.1_cds_XP_019252277.1_2282__49451"], # minor 2ODD cluster
        t["Kaladp0011s0492.1.v1.1__63787"],
    ], 
    [
        t['lcl_OX459121.1_cds_CAI9101307.1_14231__43536'], # 2ODD34
        t['rna-gnl_WGS-JAMLDZ_EVM0034906.1__4058'], 
    ], 
    [
        t['lcl_NC_039901.1_cds_XP_027112943.1_15453__13443'],
    ], 
    [
        t["rna-gnl_WGS-JAMLDZ_EVM0007978.1__4058"] # minor 2ODD cluster
    ], 
    [
        t['lcl_NC_080155.1_cds_XP_057473373.1_43266__165200'], # 2ODD35
        t['Hma1.2p1_1763F.1_g301360.1__23110'],
        t['pveT_jg21836.t1__170927'],
        t['Lj5A64G50.1__105884'],
        t['lcl_CM028326.1.g5350.t1__4392'],
        t['rna-gnl_WGS-JAMLDZ_EVM0017794.1__4058'],
        t['lcl_NC_039905.1_cds_XP_027126787.1_26673__13443'],
        t['lcl_OX459124.1_cds_CAI9114063.1_26987__43536'],
        t['lcl_PKPP01007341.1_cds_PWA53652.1_46097__35608'],
        t['rna-gnl_WGS-JAMLDZ_EVM0014154.1__4058'],
        t['lcl_OX459124.1_cds_CAI9114065.1_26989__43536'],
        t['Soltu.DM.12G010750.2__4113'],
        t['lcl_OU503036.1_cds_CAI9753676.1_1107__56036'],
        t['lcl_NC_062966.1_cds_XP_047959520.1_10496__49212'],
        t['lcl_NW_026137468.1_cds_XP_051115997.1_3322__175694']
    ], 
    [
        t['Cc02_g31130.1__49390'],  # minor 2ODD cluster
        t['Cc00_g28750.1__49390'],
        t['rna-gnl_WGS-JAMLDZ_EVM0024530.1__4058'],
        t['Cc03_g08070.1__49390'], 
    ], 
    [
        t['Anaoc.0007s0766.1__171929'],
        t['Atru_chr1_3170__47965'],
        t['FSB015489301__28930'],
        t['lcl_CM059866.1_cds_KAK0596815.1_1321__4024'],
        t['FSB011836701__28930'],
        t['CiLak.01G168500.1__32201'],
        t['Qurub.10G026200.1__3512'],
        t['Qurub.10G071600.1__3512'],
        t['lcl_CM025853.1_cds_KAB1214321.1_15709__262757'],
        t['CiLak.01G168600.1__32201'],
        t['lcl_NC_044913.1_cds_XP_030940120.1_44897__97700'],
        t['Prupe.5G048100.1__3760'],
        t['Umino00085.1__262084'],
        t['lcl_NC_024133.1_cds_XP_008243039.1_26638__102107'],
    ], 
    [
        t['lcl_NC_084796.1_cds_XP_062119410.1_20878__3486'],
        t['Umino11360.1__262084'],
        t['FvH4_3g36530.t1__57918'],
        t['lcl_NC_083380.1_cds_XP_015881543.2_5733__326968'],
        t['Kaladp0007s0029.1.v1.1__63787'],
        t['Blora12794.1__200023']
    ], 
    [
        t['VIT_209s0002g05340.1__29760'], # 2ODD36
        t['lcl_NW_017353139.1_cds_XP_018486826.1_3787__3726'],
        t['lcl_NC_045129.1_cds_XP_031385398.1_13438__22663'],
        t['lcl_JAGHRR010000128.1_cds_KAI3418995.1_3930__120290'],
        t['lcl_NW_026137602.1_cds_XP_051147854.1_28335__175694'],
        t['Hma1.2p1_0515F.1_g170445.1__23110'],
        t['lcl_CM028334.1.g444.t1__4392'],
        t['lcl_OU503038.1_cds_CAI9757760.1_5191__56036'],
        t['Vadar_g3125.t1__229202'],
        t['Acc17837.1__3625'],
        t['pveT_jg21833.t1__170927'],
        t['Lj5A64T54.1__105884'],
        t['lcl_NW_026775581.1_cds_XP_059644652.1_44213__4283'],
        t['Cc03_g10170.1__49390'],
        t['lcl_CM061492.1_cds_KAK1416820.1_25937__13708'],
        t['DCAR_016236__4039'],
        t['lcl_NC_062967.1_cds_XP_047968971.1_19558__49212'],
        t['TRINITY_DN176594_c2_g1_i2__4102'],
        t['lcl_NC_045131.1_cds_XP_031399179.1_24312__22663'],
        t['Eucgr.F02562.1__71139'],
        t['lcl_NW_006262075.1_cds_XP_006427129.2_25983__85681'],
        t['lcl_CM059866.1_cds_KAK0598226.1_1319__4024'],
        t['lcl_NC_052253.1_cds_XP_038693680.1_44458__458696'],
        t['lcl_NW_011501060.1_cds_XP_011013043.1_5386__75702'],
        t['lcl_NC_045129.1_cds_XP_031385404.1_13435__22663'],
        t['Eucgr.F02566.1__71139'],
        t['Bol038153__3712'],
        t['Lj1g0016927.1__34305'],
        t['Ler.1DRT.6g050030.1__41257'],
        t['Medtr7g090520.1__3880'],
        t['Vradi08g02800.1__157791'],
        t['Lj1g0005382.1__34305'],
        t['Vradi08g02780.1__157791'],
        t['lcl_CM059866.1_cds_KAK0596259.1_1320__4024'],
        t['Ptrif.0007s1335.1__37690'],
        t['Anaoc.0007s0767.1__171929'],
        t['Gorai.002G057700.1_PACid-26792653__29730'],
        t['Ptrif.0007s1331.1__37690'],
        t['lcl_NC_052243.1_cds_XP_038718244.1_22960__458696'],
        t['lcl_NC_083379.1_cds_XP_015868729.2_4751__326968'],
        t['lcl_NC_037092.1_cds_XP_024158425.1_33959__74649'],
        t['Umino20160.1__262084'],
        t['Umino20144.1__262084'],
        t['lcl_NC_083601.1_cds_XP_030488467.2_4391__3483'],
        t['Pparo25676.1__386216'],
        t['Dglom27811.1__34297'],
        t['Blora12130.2__200023'],
        t['lcl_NC_044908.1_cds_XP_030969963.1_22137__97700'],
        t['lcl_NC_052253.1_cds_XP_038693679.1_44455__458696'],
        t['Lus10021025_PACid-23182311__4006'],
    ], 
    [
        t['lcl_CM058159.1_cds_KAJ9692491.1_17293__103349'],
        t['Vradi08g02810.1__157791'],
        t['Ler.1DRT.1g027930.1__41257'],
        t['Lj5g0008310.1__34305'],
        t['TRINITY_DN114504_c0_g1_i1__4102'],
        t['Eucgr.J03023.1__71139'],
        t['lcl_NC_045131.1_cds_XP_031396688.1_21996__22663'],
        t['Gorai.011G079100.1_PACid-26810489__29730'],
        t['Thhalv10011669m__72664'],
        t['Anaoc.0007s0768.1__171929'],
        t['lcl_NW_006262075.1_cds_XP_006427132.2_25985__85681'],
        t['Atru_chr1_3180__47965'],
        t['Atru_chr1_3175__47965'],
        t['Potri.001G176000.1__3694'],
        t['lcl_NC_052244.1_cds_XP_038721828.1_25425__458696'],
        t['lcl_NC_065574.1_cds_XP_050234356.1_20365__3986'],
        t['lcl_NC_044908.1_cds_XP_030968259.1_22138__97700'],
        t['lcl_CM025851.1_cds_KAB1218619.1_7986__262757'],
        t['CiLak.01G242800.1__32201'],
        t['lcl_NC_037089.1_cds_XP_024183820.1_11602__74649'],
        t['lcl_NC_083379.1_cds_XP_048325586.2_4752__326968'],
        t['lcl_NC_084800.1_cds_XP_062086877.1_35687__3486'],
        t['Umino20121.1__262084'],
        t['Dglom27810.1__34297'],
        t['Tcord01333.1__703396'],
        t['Blora12127.1__200023'],
    ], 
    [
        t['CiLak.02G150500.1__32201'],
        t['FSB014894201__28930'],
        t['lcl_NC_044908.1_cds_XP_030968569.1_22136__97700']
 ], 
 [
        t['lcl_PKPP01002761.1_cds_PWA73280.1_26522__35608'], # 2ODD35
        t['rna-gnl_WGS-JAMLDZ_EVM0012381.1__4058'],
        t['lcl_CM028334.1.g446.t1__4392'],
        t['Lj5A64T53.1__105884'],
        t['Hma1.2p1_3098F.1_g332680.1__23110'],
        t['lcl_NW_026775581.1_cds_XP_059639306.1_44215__4283'],
        t['Vadar_g3126.t1__229202'],
        t['Acc12585.1__3625'],
        t['pveT_jg21835.t1__170927'],
        t['Pg_S0988.1__4054'],
        t['DCAR_016235__4039'],
        t['lcl_OX459124.1_cds_CAI9114062.1_26986__43536'],
        t['lcl_OU503036.1_cds_CAI9753675.1_1106__56036'],
        t['lcl_NC_062966.1_cds_XP_047959442.1_10495__49212'],
        t['lcl_NW_026137602.1_cds_XP_051147853.1_28334__175694'],
        t['DCAR_007548__4039'],
        t['lcl_CM061487.1_cds_KAK1434086.1_11004__13708'],
        t['HanXRQChr09g0275101__4232'],
        t['Lj5C65T2.1__105884'],
        t['lcl_NW_026775581.1_cds_XP_059644765.1_44212__4283'],
        t['lcl_NC_080143.1_cds_XP_057503960.1_20127__165200'],
        t['pveT_jg21832.t1__170927'],
        t['lcl_CM028334.1.g443.t1__4392'],
        t['rna-gnl_WGS-JAMLDZ_EVM0018858.1__4058'],
        t['SMEL_006g260960.1.01__4111'],
        t['lcl_NC_039905.1_cds_XP_027126923.1_26675__13443'],
        t['lcl_OU503036.1_cds_CAI9753673.1_1104__56036'],
        t['MSTRG.22708.1__82927'],
        t['lcl_NC_062966.1_cds_XP_047960078.1_10492__49212'],
        t['lcl_NW_026137587.1_cds_XP_051135949.1_16701__175694'],
        t['lcl_NW_026137587.1_cds_XP_051135915.1_16702__175694']
    ], 
    [
        t['lcl_CM008454.1_cds_PHT32940.1_29277__33114'], # 2ODD38
        t['lcl_CM008444.1_cds_PHT55775.1_4261__33114'],
        t['TRINITY_DN232929_c4_g6_i3__4102'],
        t['lcl_CM008445.1_cds_PHT51893.1_6355__33114'],
        t['lcl_NC_083337.1_cds_XP_060201071.1_373__112863'],
        t['lcl_NC_083338.1_cds_XP_060200949.1_5353__112863'],
        t['Soltu.DM.01G002190.1__4113'],
        t['SMEL_001g152420.1.01__4111'],
        t['lcl_NC_028638.1_cds_XP_015067085.1_5485__28526'],
        t['TRINITY_DN725607_c0_g1_i1__189803'],
        t['Vadar_g43141.t1__229202'],
        t['KT390173__DPS__etoposide_biosynthesis__93608'],
        t['lcl_NW_026775571.1_cds_XP_059664661.1_3040__4283'],
        t['AF417859__AOP3__glucosinolate_biosynthesis__3702'],
        t['AF417858__AOP2__glucosinolate_biosynthesis__3702'],
        t['lcl_CM042048.1_cds_KAI3759661.1_7640__4217'],
        t['lcl_NW_026137594.1_cds_XP_051141059.1_21911__175694'],
        t['lcl_NC_083379.1_cds_XP_048320867.2_1719__326968'],
        t['Umino20164.1__262084'],
        t['Cc08_g04750.1__49390'],
        t['Lus10023024_PACid-23146145__4006'],
        t['lcl_NC_045131.1_cds_XP_031396066.1_23250__22663'],
        t['lcl_NC_052241.1_cds_XP_038713292.1_19041__458696'],
        t['Gobar.A11G255200.1__3634'],
        t['Medtr4g011690.1__3880'],
        t['Ca_13466__3827'],
        t['Ler.1DRT.4g003660.1__41257'],
        t['Ca_13467__3827'],
        t['Lj3g0013768.1__34305'],
        t['lcl_OX451740.1_cds_CAI8613451.1_28042__3906'],
        t['Lj1g0004205.1__34305'],
        t['Vradi0154s00060.1__157791'],
        t['lcl_NC_052245.1_cds_XP_038722105.1_29001__458696'],
        t['lcl_NC_065574.1_cds_XP_050231702.1_19581__3986'],
        t['Lus10005516_PACid-23158832__4006'],
        t['lcl_NC_052393.1_cds_XP_008784635.1_5049__42345'],
        t['Vadar_g3127.t1__229202'],
        t['VIT_209s0002g05280.1__29760'],
        t['Acc17838.1__3625'],
        t['lcl_NW_026775581.1_cds_XP_059643089.1_44218__4283'],
        t['lcl_CM028334.1.g447.t1__4392'],
        t['Lj9C505T3.1__105884'],
        t['Pg_S0213.71__4054'],
        t['DCAR_016234__4039'],
        t['lcl_CM027379.1_cds_KAF9612984.1_11656__261450'],
        t['lcl_PKPP01007341.1_cds_PWA53651.1_46096__35608'],
        t['Hma1.2p1_0515F.1_g170455.1__23110'],
        t['Ah1G51780.1__81970'],
        t['Cvari18296.1__869952'],
        t['Umino20166.1__262084'],
        t['lcl_NC_083379.1_cds_XP_015868730.2_4749__326968'],
        t['lcl_NC_044908.1_cds_XP_030973050.1_22133__97700'],
        t['CiLak.01G242700.1__32201'],
        t['lcl_CM025851.1_cds_KAB1218616.1_7983__262757'],
        t['Prupe.3G075700.1__3760'],
        t['Dglom27813.1__34297'],
        t['Blora09835.1__200023'],
        t['lcl_NW_006262075.1_cds_XP_006427127.2_25981__85681'],
        t['Potri.001G176200.1__3694'],
        t['lcl_NC_065574.1_cds_XP_050231462.1_20561__3986'],
        t['lcl_JAGHRR010000039.1_cds_KAI3437248.1_15716__120290'],
        t['Anaoc.0007s0764.1__171929'],
        t['lcl_CM059866.1_cds_KAK0595174.1_1315__4024'],
        t['OMO54201__210143'],
        t['lcl_NC_052253.1_cds_XP_038694945.1_44457__458696'],
        t['Lus10016659_PACid-23143917__4006'],
        t['pveT_jg11526.t1__170927'],
        t['lcl_NC_083601.1_cds_XP_060957775.1_4387__3483'],
        t['lcl_OX459123.1_cds_CAI9108625.1_21549__43536'],
        t['Kaladp0008s0237.1.v1.1__63787'],
        t['Soltu.DM.06G023510.1__4113'],
        t['rna-gnl_WGS-JAMLDZ_EVM0014967.1__4058'],
        t['lcl_OU503036.1_cds_CAI9753677.1_1108__56036'],
        t['MSTRG.11101.1__82927'],
        t['lcl_NC_062966.1_cds_XP_047958340.1_10499__49212'],
        t['lcl_NC_062965.1_cds_XP_047956864.1_861__49212'], 
    ], 
    [
        t['Mba05_g08120.1.v1.1__52838'],
    ], 
    [
        t['CALSI_Maker00003989__746888'],
        t['Mba03_g21520.1.v1.1__52838'],
        t['evm.model.AsparagusV1_08.662__4686'],
        t['Ma04_t25080.1__4641'],
        t['lcl_CP136891.1_cds_WOK98631.1_7322__4628'],
        t['Mba05_g08100.1.v1.1__52838'],
        t['Mtr04_g32474.1__320322'],
        t['lcl_CP136893.1_cds_WOL03332.1_12022__4628'],
        t['CALSI_Maker00009131__746888'],
        t['ORUFI05G28750.1__4529'],
        t['ORUFI11G09650.1__4529'],
    ], 
    [
        t['lcl_NC_052396.1_cds_XP_008794218.1_12468__42345'],
        t['Spipo4G0051700__29656'],
    ], 
    [
        t['lcl_CP136890.1_cds_WOK92059.1_750__4628'],
        t['CALSI_Maker00050580__746888'],
        t['evm.model.AsparagusV1_01.591__4686'],
        t['ASH_rna4294__1088818'],
    ], 
    [
        t['strangu_020525-RA__38733']
    ], 
    [
        t['lcl_NC_053024.1_cds_XP_048567703.1_14633__4572']  # minor 2ODD cluster
    ]
]

def test_get_clusters_without_candidates():
    result = get_clusters(t)
    assert result == expected_clusters_without_candidates


# take the test tree with all input sequences, manually select candidates, and create a dictionary with manually assigned candidates as keys that map to the expected output of the annotation service.
# map candidate name to expected cluster id resulted from the annotation service.
# The expected cluster id IS NOT the cluster id of the actual input sequence 
expected_clusters_with_candidates = [
    [
        t["lcl_NC_084852.1_cds_XP_061960601.1_3168__3691"] # candidate
    ], 
    [
        t["rna-XM_031626300.2__210225"], # 2ODD34
    ], 
    [
        t["CKAN_00595100__337451"], 
        t["Aristolochia_fimbriatalcl_CM034078.1_cds_KAG9459387.1_3895__158543"], 
        t["Gobar.D13G239800.1__3634"],
        t["Lj8A611G43.1__105884"], 
        t["FvH4_2g29990.t1__57918"],
        t['Casgl23S04257__3522'],
        t['CiLak.07G220400.1__32201'],
        t['lcl_NC_065570.1_cds_XP_050209022.1_2785__3986'],
        t['Lus10013130_PACid-23164132__4006'],
    ], 
    [
        t['evm_27.model.AmTr_v1.0_scaffold00174.4__13333'],
    ], 
    [
        t['RZC71586__3469'],
    ], 
    [
        t['lcl_NW_010729080.1_cds_XP_010249917.1_5454__4432'],
        t['rna-XM_042628350.1__60698'],
        t['lcl_CM056811.1_cds_KAJ8637692.1_11945__3435'],
        t['Aristolochia_fimbriatalcl_CM034078.1_cds_KAG9459386.1_3894__158543'],
    ], 
    [
        t['Kaladp1006s0012.1.v1.1__63787'],
        t['VIT_204s0008g04920.2__29760'], # candidate
        t['Ptrif.0001s0318.1__37690'],# candidate
        t['Atru_chr7_2342__47965'],# candidate
        t['lcl_NW_019168159.1_cds_XP_022728803.1_47147__66656'],# candidate
        t['lcl_NC_045134.1_cds_XP_031407529.1_33286__22663'],# candidate
        t['GWHTACBH000860__413952'],# candidate
        t['Potri.006G248000.1__3694'],# candidate
        t['lcl_NC_065570.1_cds_XP_050209042.1_2787__3986'],# candidate
        t['Lus10008097_PACid-23169790__4006'],# candidate
        t['Acc06446.1__3625'],
        t['Acc16585.1__3625'],
        t['pveT_jg24143.t1__170927'],
        t['pveT_jg28701.t1__170927'],
        t['lcl_NC_083383.1_cds_XP_048332159.2_17459__326968'],
        t['Cnepa31456.1__79760'],
        t['Umino40998.1__262084'],
        t['FvH4_2g29991.t1__57918'],
        t['Lj1g0022518.1__34305'],
        t['Lj4g0022664.1__34305'],
        t['Medtr1g011600.1__3880'],
        t['Ca_10933__3827'],
        t['Ler.1DRT.1g009880.1__41257'],
        t['Ler.1DRT.3g073930.1__41257'],
        t['Medtr3g108520.1__3880'],
        t['Ca_23071__3827'],
        t['Ca_01791__3827'],
        t['Medtr4g021360.1__3880'],
        t['Ca_05192__3827'],
        t['Medtr4g021380.1__3880'],
        t['lcl_OX451737.1_cds_CAI8600332.1_14923__3906']
    ], 
    [
        t["pveT_jg20218.t1__170927"], # minor 2ODD cluster
        t["pveT_jg20220.t1__170927"],
        t["Acc09929.1__3625"], # candidate
        t["pveT_jg33784.t1__170927"] # minor 2ODD cluster
    ], 
    [
        t['Eucgr.C03724.1__71139'], # 2ODD38
        t['lcl_NC_045130.1_cds_XP_031390300.1_19564__22663'],
        t['Eucgr.C04147.1__71139'],
    ], 
    [
        t['Acc23206.1__3625'],
        t['lcl_CM027379.1_cds_KAF9612985.1_11655__261450'],
        t['RZC71605__3469'],
        t['rna-XM_042647017.1__60698'],
    ], 
    [
        t['lcl_CM027379.1_cds_KAF9612671.1_9820__261450']
    ],
    [
        t['Lj7A541T77.1__105884'] # minor 2ODD cluster
    ],
    [
        t['Pg_S7309.1__4054'],  # 2ODD34
        t['DCAR_015859__4039'],
    ], 
    [
        t['CsatW809539.1__3659'], 
    ], 
    [
        t['lcl_NW_026775580.1_cds_XP_059633894.1_39517__4283'],
        t['Dglom27812.1__34297'],
        t['lcl_CM059866.1_cds_KAK0595785.1_1317__4024'],
        t['Ptrif.0007s1329.1__37690'],
        t['Potri.001G176500.1__3694'],
        t['lcl_NC_065574.1_cds_XP_050234707.1_20560__3986'],
        t['Umino20162.1__262084'],
        t['lcl_NC_083379.1_cds_XP_015868728.2_4750__326968'],
        t['lcl_NC_084800.1_cds_XP_062086870.1_35676__3486'],
        t['lcl_NC_045131.1_cds_XP_031396687.1_21997__22663'],
        t['Eucgr.J03022.1__71139'],
        t['VIT_209s0002g05290.1__29760'],
        t['lcl_NW_026775580.1_cds_XP_059636056.1_39514__4283'],
        t['Acc08999.1__3625'],
        t['Acc07679.1__3625'],
        t['FSB010835101__28930'],
        t['DCAR_015858__4039']
    ], 
    [
        t['Casgl236S24615__3522'], # minor 2ODD cluster
        t['CiLak.15G093300.1__32201'],
        t['Qurub.10G093800.1__3512'],
        t['lcl_CM025849.1_cds_KAB1226463.1_1752__262757'],
        t['Tcord05348.2__703396'],
        t['HU01G00826.1__176265'],
        t['Sc13g0006470.01__3999'],
        t['TRINITY_DN108901_c0_g1_i1__223224'],
        t['CM008280.1.CM008280.1.g987.t1__62330']
    ], 
    [
        t['pveT_jg9723.t1__170927'] # 2ODD34
    ],
    [
        t["Cc06_g15940.1__49390"], # minor 2ODD cluster
    ], 
    [
        t["Gobar.D05G137100.1__3634"],
        t["lcl_OX459118.1_cds_CAI9089985.1_2909__43536"],
    ],
    [
        t['TRINITY_DN197356_c3_g1_i1__4102'], # 2ODD34
        t['lcl_CM061488.1_cds_KAK1430368.1_13046__13708'],
        t['Oeu000043.1__4146'],
        t['lcl_NW_026137546.1_cds_XP_051123007.1_5823__175694'],
        t['lcl_NW_025952766.1_cds_XP_047954044.1_46391__49212'],
        t['lcl_NC_028637.1_cds_XP_015064474.1_2338__28526'],
        t['Cc00_g25210.1__49390'],
        t['rna-gnl_WGS-JAMLDZ_EVM0036354.1__4058'],
        t['lcl_OX459121.1_cds_CAI9101309.1_14233__43536']
    ],
    [
        t["lcl_NC_031989.1_cds_XP_019252277.1_2282__49451"],# candidate
        t["Kaladp0011s0492.1.v1.1__63787"],# candidate
        t['lcl_OX459121.1_cds_CAI9101307.1_14231__43536'], # candidate
        t['rna-gnl_WGS-JAMLDZ_EVM0034906.1__4058'], # candidate
        t['lcl_NC_039901.1_cds_XP_027112943.1_15453__13443'],# candidate
        t["rna-gnl_WGS-JAMLDZ_EVM0007978.1__4058"],# candidate
        t['lcl_NC_080155.1_cds_XP_057473373.1_43266__165200'], # 2ODD35
        t['Hma1.2p1_1763F.1_g301360.1__23110'],
        t['pveT_jg21836.t1__170927'],
        t['Lj5A64G50.1__105884'],
        t['lcl_CM028326.1.g5350.t1__4392'],
        t['rna-gnl_WGS-JAMLDZ_EVM0017794.1__4058'],
        t['lcl_NC_039905.1_cds_XP_027126787.1_26673__13443'],
        t['lcl_OX459124.1_cds_CAI9114063.1_26987__43536'],
        t['lcl_PKPP01007341.1_cds_PWA53652.1_46097__35608'],
        t['rna-gnl_WGS-JAMLDZ_EVM0014154.1__4058'],
        t['lcl_OX459124.1_cds_CAI9114065.1_26989__43536'],
        t['Soltu.DM.12G010750.2__4113'],
        t['lcl_OU503036.1_cds_CAI9753676.1_1107__56036'],
        t['lcl_NC_062966.1_cds_XP_047959520.1_10496__49212'],
        t['lcl_NW_026137468.1_cds_XP_051115997.1_3322__175694'],
    ], 
    [
        t['Cc02_g31130.1__49390'],  # minor 2ODD cluster
        t['Cc00_g28750.1__49390'],
        t['rna-gnl_WGS-JAMLDZ_EVM0024530.1__4058'],
        t['Cc03_g08070.1__49390'], 
    ], 
    [
        t['Anaoc.0007s0766.1__171929'],
        t['Atru_chr1_3170__47965'],
        t['FSB015489301__28930'],
        t['lcl_CM059866.1_cds_KAK0596815.1_1321__4024'],
        t['FSB011836701__28930'],
        t['CiLak.01G168500.1__32201'],
        t['Qurub.10G026200.1__3512'],
        t['Qurub.10G071600.1__3512'],
        t['lcl_CM025853.1_cds_KAB1214321.1_15709__262757'],
        t['CiLak.01G168600.1__32201'],
        t['lcl_NC_044913.1_cds_XP_030940120.1_44897__97700'],
        t['Prupe.5G048100.1__3760'],
        t['Umino00085.1__262084'],
        t['lcl_NC_024133.1_cds_XP_008243039.1_26638__102107'],
    ], 
    [
        t['lcl_NC_084796.1_cds_XP_062119410.1_20878__3486'],
        t['Umino11360.1__262084'],
        t['FvH4_3g36530.t1__57918'],
        t['lcl_NC_083380.1_cds_XP_015881543.2_5733__326968'],
        t['Kaladp0007s0029.1.v1.1__63787'],
        t['Blora12794.1__200023']
    ], 
    [
        t['VIT_209s0002g05340.1__29760'], # 2ODD36
        t['lcl_NW_017353139.1_cds_XP_018486826.1_3787__3726'],
        t['lcl_NC_045129.1_cds_XP_031385398.1_13438__22663'],
        t['lcl_JAGHRR010000128.1_cds_KAI3418995.1_3930__120290'],
        t['lcl_NW_026137602.1_cds_XP_051147854.1_28335__175694'],
        t['Hma1.2p1_0515F.1_g170445.1__23110'],
        t['lcl_CM028334.1.g444.t1__4392'],
        t['lcl_OU503038.1_cds_CAI9757760.1_5191__56036'],
        t['Vadar_g3125.t1__229202'],
        t['Acc17837.1__3625'],
        t['pveT_jg21833.t1__170927'],
        t['Lj5A64T54.1__105884'],
        t['lcl_NW_026775581.1_cds_XP_059644652.1_44213__4283'],
        t['Cc03_g10170.1__49390'],
        t['lcl_CM061492.1_cds_KAK1416820.1_25937__13708'],
        t['DCAR_016236__4039'],
        t['lcl_NC_062967.1_cds_XP_047968971.1_19558__49212'],
        t['TRINITY_DN176594_c2_g1_i2__4102'],
        t['lcl_NC_045131.1_cds_XP_031399179.1_24312__22663'],
        t['Eucgr.F02562.1__71139'],
        t['lcl_NW_006262075.1_cds_XP_006427129.2_25983__85681'],
        t['lcl_CM059866.1_cds_KAK0598226.1_1319__4024'],
        t['lcl_NC_052253.1_cds_XP_038693680.1_44458__458696'],
        t['lcl_NW_011501060.1_cds_XP_011013043.1_5386__75702'],
        t['lcl_NC_045129.1_cds_XP_031385404.1_13435__22663'],
        t['Eucgr.F02566.1__71139'],
        t['Bol038153__3712'],   # candidate
        t['Lj1g0016927.1__34305'],
        t['Ler.1DRT.6g050030.1__41257'],
        t['Medtr7g090520.1__3880'],
        t['Vradi08g02800.1__157791'],
        t['Lj1g0005382.1__34305'],
        t['Vradi08g02780.1__157791'],
        t['lcl_CM059866.1_cds_KAK0596259.1_1320__4024'],
        t['Ptrif.0007s1335.1__37690'],
        t['Anaoc.0007s0767.1__171929'],
        t['Gorai.002G057700.1_PACid-26792653__29730'],
        t['Ptrif.0007s1331.1__37690'],
        t['lcl_NC_052243.1_cds_XP_038718244.1_22960__458696'],
        t['lcl_NC_083379.1_cds_XP_015868729.2_4751__326968'],
        t['lcl_NC_037092.1_cds_XP_024158425.1_33959__74649'],
        t['Umino20160.1__262084'],
        t['Umino20144.1__262084'],
        t['lcl_NC_083601.1_cds_XP_030488467.2_4391__3483'],
        t['Pparo25676.1__386216'],
        t['Dglom27811.1__34297'],
        t['Blora12130.2__200023'],
        t['lcl_NC_044908.1_cds_XP_030969963.1_22137__97700'],
        t['lcl_NC_052253.1_cds_XP_038693679.1_44455__458696'],
        t['Lus10021025_PACid-23182311__4006'],
    ], 
    [
        t['lcl_CM058159.1_cds_KAJ9692491.1_17293__103349'],
        t['Vradi08g02810.1__157791'],
        t['Ler.1DRT.1g027930.1__41257'],
        t['Lj5g0008310.1__34305'],
        t['TRINITY_DN114504_c0_g1_i1__4102'],
        t['Eucgr.J03023.1__71139'],
        t['lcl_NC_045131.1_cds_XP_031396688.1_21996__22663'],
        t['Gorai.011G079100.1_PACid-26810489__29730'],
        t['Thhalv10011669m__72664'],
        t['Anaoc.0007s0768.1__171929'],
        t['lcl_NW_006262075.1_cds_XP_006427132.2_25985__85681'],
        t['Atru_chr1_3180__47965'],
        t['Atru_chr1_3175__47965'],
        t['Potri.001G176000.1__3694'],
        t['lcl_NC_052244.1_cds_XP_038721828.1_25425__458696'],
        t['lcl_NC_065574.1_cds_XP_050234356.1_20365__3986'],
        t['lcl_NC_044908.1_cds_XP_030968259.1_22138__97700'],
        t['lcl_CM025851.1_cds_KAB1218619.1_7986__262757'],
        t['CiLak.01G242800.1__32201'],
        t['lcl_NC_037089.1_cds_XP_024183820.1_11602__74649'],
        t['lcl_NC_083379.1_cds_XP_048325586.2_4752__326968'],
        t['lcl_NC_084800.1_cds_XP_062086877.1_35687__3486'],
        t['Umino20121.1__262084'],
        t['Dglom27810.1__34297'],
        t['Tcord01333.1__703396'],
        t['Blora12127.1__200023'],
    ], 
    [
        t['CiLak.02G150500.1__32201'],
        t['FSB014894201__28930'],
        t['lcl_NC_044908.1_cds_XP_030968569.1_22136__97700']
    ], 
    [
        t['lcl_PKPP01002761.1_cds_PWA73280.1_26522__35608'], # 2ODD35
        t['rna-gnl_WGS-JAMLDZ_EVM0012381.1__4058'],
        t['lcl_CM028334.1.g446.t1__4392'],
        t['Lj5A64T53.1__105884'],
        t['Hma1.2p1_3098F.1_g332680.1__23110'],
        t['lcl_NW_026775581.1_cds_XP_059639306.1_44215__4283'],
        t['Vadar_g3126.t1__229202'],
        t['Acc12585.1__3625'],
        t['pveT_jg21835.t1__170927'],
        t['Pg_S0988.1__4054'],
        t['DCAR_016235__4039'],
        t['lcl_OX459124.1_cds_CAI9114062.1_26986__43536'],
        t['lcl_OU503036.1_cds_CAI9753675.1_1106__56036'],
        t['lcl_NC_062966.1_cds_XP_047959442.1_10495__49212'],
        t['lcl_NW_026137602.1_cds_XP_051147853.1_28334__175694'],
        t['DCAR_007548__4039'],
        t['lcl_CM061487.1_cds_KAK1434086.1_11004__13708'],
        t['HanXRQChr09g0275101__4232'],
        t['Lj5C65T2.1__105884'],
        t['lcl_NW_026775581.1_cds_XP_059644765.1_44212__4283'],
        t['lcl_NC_080143.1_cds_XP_057503960.1_20127__165200'],
        t['pveT_jg21832.t1__170927'],
        t['lcl_CM028334.1.g443.t1__4392'],
        t['rna-gnl_WGS-JAMLDZ_EVM0018858.1__4058'],
        t['SMEL_006g260960.1.01__4111'],
        t['lcl_NC_039905.1_cds_XP_027126923.1_26675__13443'],
        t['lcl_OU503036.1_cds_CAI9753673.1_1104__56036'],
        t['MSTRG.22708.1__82927'],
        t['lcl_NC_062966.1_cds_XP_047960078.1_10492__49212'],
        t['lcl_NW_026137587.1_cds_XP_051135949.1_16701__175694'],
        t['lcl_NW_026137587.1_cds_XP_051135915.1_16702__175694'],
    ], 
    [
        t['lcl_CM008454.1_cds_PHT32940.1_29277__33114'], # 2ODD38
        t['lcl_CM008444.1_cds_PHT55775.1_4261__33114'],
        t['TRINITY_DN232929_c4_g6_i3__4102'],
        t['lcl_CM008445.1_cds_PHT51893.1_6355__33114'],
        t['lcl_NC_083337.1_cds_XP_060201071.1_373__112863'],
        t['lcl_NC_083338.1_cds_XP_060200949.1_5353__112863'],
        t['Soltu.DM.01G002190.1__4113'],
        t['SMEL_001g152420.1.01__4111'],
        t['lcl_NC_028638.1_cds_XP_015067085.1_5485__28526'],
        t['TRINITY_DN725607_c0_g1_i1__189803'],
        t['Vadar_g43141.t1__229202'],
        t['KT390173__DPS__etoposide_biosynthesis__93608'],
        t['lcl_NW_026775571.1_cds_XP_059664661.1_3040__4283'],
        t['AF417859__AOP3__glucosinolate_biosynthesis__3702'],
        t['AF417858__AOP2__glucosinolate_biosynthesis__3702'],
        t['lcl_CM042048.1_cds_KAI3759661.1_7640__4217'],
        t['lcl_NW_026137594.1_cds_XP_051141059.1_21911__175694'],
        t['lcl_NC_083379.1_cds_XP_048320867.2_1719__326968'],
        t['Umino20164.1__262084'],
        t['Cc08_g04750.1__49390'],
        t['Lus10023024_PACid-23146145__4006'],
        t['lcl_NC_045131.1_cds_XP_031396066.1_23250__22663'],
        t['lcl_NC_052241.1_cds_XP_038713292.1_19041__458696'],
        t['Gobar.A11G255200.1__3634'],
        t['Medtr4g011690.1__3880'],
        t['Ca_13466__3827'],
        t['Ler.1DRT.4g003660.1__41257'],
        t['Ca_13467__3827'],
        t['Lj3g0013768.1__34305'],
        t['lcl_OX451740.1_cds_CAI8613451.1_28042__3906'],
        t['Lj1g0004205.1__34305'],
        t['Vradi0154s00060.1__157791'],
        t['lcl_NC_052245.1_cds_XP_038722105.1_29001__458696'],
        t['lcl_NC_065574.1_cds_XP_050231702.1_19581__3986'],
        t['Lus10005516_PACid-23158832__4006'],
        t['lcl_NC_052393.1_cds_XP_008784635.1_5049__42345'],
        t['Vadar_g3127.t1__229202'],
        t['VIT_209s0002g05280.1__29760'],
        t['Acc17838.1__3625'],
        t['lcl_NW_026775581.1_cds_XP_059643089.1_44218__4283'],
        t['lcl_CM028334.1.g447.t1__4392'],
        t['Lj9C505T3.1__105884'],
        t['Pg_S0213.71__4054'],
        t['DCAR_016234__4039'],
        t['lcl_CM027379.1_cds_KAF9612984.1_11656__261450'],
        t['lcl_PKPP01007341.1_cds_PWA53651.1_46096__35608'],
        t['Hma1.2p1_0515F.1_g170455.1__23110'],
        t['Ah1G51780.1__81970'],
        t['Cvari18296.1__869952'],
        t['Umino20166.1__262084'],
        t['lcl_NC_083379.1_cds_XP_015868730.2_4749__326968'],
        t['lcl_NC_044908.1_cds_XP_030973050.1_22133__97700'],
        t['CiLak.01G242700.1__32201'],
        t['lcl_CM025851.1_cds_KAB1218616.1_7983__262757'],
        t['Prupe.3G075700.1__3760'],
        t['Dglom27813.1__34297'],
        t['Blora09835.1__200023'],
        t['lcl_NW_006262075.1_cds_XP_006427127.2_25981__85681'],
        t['Potri.001G176200.1__3694'],
        t['lcl_NC_065574.1_cds_XP_050231462.1_20561__3986'],
        t['lcl_JAGHRR010000039.1_cds_KAI3437248.1_15716__120290'],
        t['Anaoc.0007s0764.1__171929'],
        t['lcl_CM059866.1_cds_KAK0595174.1_1315__4024'],
        t['OMO54201__210143'],
        t['lcl_NC_052253.1_cds_XP_038694945.1_44457__458696'],
        t['Lus10016659_PACid-23143917__4006'],
        t['pveT_jg11526.t1__170927'],
        t['lcl_NC_083601.1_cds_XP_060957775.1_4387__3483'],
        t['lcl_OX459123.1_cds_CAI9108625.1_21549__43536'],
        t['Kaladp0008s0237.1.v1.1__63787'],
        t['Soltu.DM.06G023510.1__4113'],
        t['rna-gnl_WGS-JAMLDZ_EVM0014967.1__4058'],
        t['lcl_OU503036.1_cds_CAI9753677.1_1108__56036'],
        t['MSTRG.11101.1__82927'],
        t['lcl_NC_062966.1_cds_XP_047958340.1_10499__49212'],
        t['lcl_NC_062965.1_cds_XP_047956864.1_861__49212'], 
        t['Mba05_g08120.1.v1.1__52838'],
        t['CALSI_Maker00003989__746888'],
        t['Mba03_g21520.1.v1.1__52838'],
        t['evm.model.AsparagusV1_08.662__4686'],
        t['Ma04_t25080.1__4641'],
        t['lcl_CP136891.1_cds_WOK98631.1_7322__4628'],
        t['Mba05_g08100.1.v1.1__52838'],
        t['Mtr04_g32474.1__320322'],
        t['lcl_CP136893.1_cds_WOL03332.1_12022__4628'],
        t['CALSI_Maker00009131__746888'],
        t['ORUFI05G28750.1__4529'],
        t['ORUFI11G09650.1__4529'],
        t['lcl_NC_052396.1_cds_XP_008794218.1_12468__42345'],
        t['Spipo4G0051700__29656'],
        t['lcl_CP136890.1_cds_WOK92059.1_750__4628'],
        t['CALSI_Maker00050580__746888'],
        t['evm.model.AsparagusV1_01.591__4686'],
        t['ASH_rna4294__1088818'],
        t['strangu_020525-RA__38733'],
        t['lcl_NC_053024.1_cds_XP_048567703.1_14633__4572'], # candidate
    ]
]

def test_clusters_with_candidates():
    manual_candidates = {
    "lcl_NC_084852.1_cds_XP_061960601.1_3168__3691" : "2ODD34",
    'VIT_204s0008g04920.2__29760' : '2ODD34',
    'Ptrif.0001s0318.1__37690' : '2ODD34', 
    'Atru_chr7_2342__47965' : '2ODD34',
    'lcl_NW_019168159.1_cds_XP_022728803.1_47147__66656' : '2ODD34',
    'lcl_NC_045134.1_cds_XP_031407529.1_33286__22663' : '2ODD34',
    'GWHTACBH000860__413952' : '2ODD34',
    'Potri.006G248000.1__3694' : '2ODD34',
    'lcl_NC_065570.1_cds_XP_050209042.1_2787__3986' : '2ODD34',
    'Lus10008097_PACid-23169790__4006' : '2ODD34', 
    'Acc09929.1__3625': "minor_2ODD_cluster",
    'lcl_NC_031989.1_cds_XP_019252277.1_2282__49451': "2ODD35",
    'Kaladp0011s0492.1.v1.1__63787': "2ODD35",
    'lcl_OX459121.1_cds_CAI9101307.1_14231__43536': "2ODD35",
    'rna-gnl_WGS-JAMLDZ_EVM0034906.1__4058': "2ODD35",
    'lcl_NC_039901.1_cds_XP_027112943.1_15453__13443': "2ODD35",
    'rna-gnl_WGS-JAMLDZ_EVM0007978.1__4058': "2ODD35", 
    'Bol038153__3712': "2ODD36",
    "lcl_NC_053024.1_cds_XP_048567703.1_14633__4572": "unresolved"
    }

    for leaf in t:
        if leaf.name in manual_candidates:
            leaf.props["two_odd_id"] = "candidate"

    result = get_clusters(t)
    assert result == expected_clusters_with_candidates



#%% ====================TEST CLUSTERING WITH NESTED 2ODD IDs====================

# test tree with nested 2ODD IDs (e.g. 2ODD14A nested within 2ODD14)
test_nested_tree_path = Path(__file__).parents[1] / "data" / "test_nested_tree.nwk"
t_nested = PhyloTree(open(test_nested_tree_path), sp_naming_function=lambda name: name.split('__')[-1])
tax2names, tax2lineages, tax2rank = t_nested.annotate_ncbi_taxa(taxid_attr='species')
candidate_headers, ingroup_headers = split_seqs_by_2ODD_membership(seq_to_2ODD_id=seq_to_2ODD_id, tree=t_nested)
assign_2ODD_props(
    tree=t_nested, seq_to_2ODD_id=seq_to_2ODD_id, 
    candidate_headers=candidate_headers
)
assign_plant_group_props(tree=t_nested)
dist_dict_nested = build_distance_lookup(t_nested)


# explore_2ODD_IDs(t_nested)
#%%

def test_get_clusters_nested_without_candidates():
    manual_nested_seqs_12B = [
        t_nested['Vradi08g06870.1__157791'],
        t_nested['Medtr1g070135.1__3880'],
        t_nested['Ler.1DRT.1g071540.1__41257'],
        t_nested['Medtr1g070120.1__3880'],
        t_nested['Ler.1DRT.1g071570.1__41257'],
        t_nested['Medtr1g070080.1__3880'],
        t_nested['Ler.1DRT.1g071590.1__41257'],
        t_nested['Medtr1g070085.3__3880'],
    ]

    manual_nested_seqs_12A = [
        t_nested['Cc01_g06970.1__49390'],
        t_nested['Cc01_g09870.1__49390'],
        t_nested['rna-gnl_WGS-JAMLDZ_EVM0004127.1__4058'],
    ]

    for leaf in t_nested:
        if leaf in manual_nested_seqs_12B:
            leaf.props["two_odd_id"] = "2ODD12B"
        elif leaf in manual_nested_seqs_12A:
            leaf.props["two_odd_id"] = "2ODD12A"

    expected_clusters_nestd_without_candidates = [

        [
            t_nested['Kaladp0089s0096.1.v1.1__63787'],
            t_nested['lcl_CM027381.1_cds_KAF9605627.1_22189__261450'],
            t_nested['Vadar_g6879.t1__229202'],
            t_nested['Acc25611.1__3625'],
            t_nested['Dglom04183.1__34297'],
            t_nested['Anaoc.0002s1042.1__171929'],
            t_nested['Eucgr.B01535.1__71139'],
            t_nested['Qurub.09G168900.1__3512'],
            t_nested['Prupe.8G225800.1__3760'],
            t_nested['Casgl267S13418__3522'],
            t_nested['Casgl528S19423__3522'],
            t_nested['lcl_CM025853.1_cds_KAB1213729.1_15117__262757'],
            t_nested['lcl_NC_083604.1_cds_XP_030496368.2_17710__3483'],
            t_nested['Lj5g0001478.1__34305'],
            t_nested['Lus10035782_PACid-23147456__4006'],
            t_nested['lcl_NW_011499926.1_cds_XP_011038934.1_22618__75702'],
            t_nested['Gobar.D08G228100.1__3634'],
            t_nested['Anaoc.1020s0004.1__171929'],
            t_nested['Ptrif.0003s3349.1__37690'],
            t_nested['lcl_CM059870.1_cds_KAK0580988.1_15865__4024'],
            t_nested['FvH4_2g21630.t4__57918'],
            t_nested['Sivu_ALN-8597__42043'],
            t_nested['Sila-35377__37657'],
            t_nested['Sc24g0003120.01__3999'],
            t_nested['TRINITY_DN229963_c3_g1_i1__122310'],
            t_nested['TRINITY_DN105553_c0_g1_i1__223224'],
            t_nested['Sc17g0000710.01__3999'],
            t_nested['DIACA2-20392__3570'],
            t_nested['Sivu_ALN-38268__42043'],
            t_nested['CM029399.1.CM029399.1.g1105.t1__3670'],
            t_nested['lcl_CM059870.1_cds_KAK0581506.1_15868__4024'],
            t_nested['lcl_NW_006262339.1_cds_XP_024044614.1_12214__85681'],
            t_nested['Thecc.01G260000.1__3641'],
            t_nested['Dglom04181.1__34297'],
            t_nested['Blora11689.1__200023'],
            t_nested['Qurub.09G177400.1__3512'],
            t_nested['lcl_NC_083387.1_cds_XP_015891362.3_28814__326968'],
            t_nested['lcl_NC_083604.1_cds_XP_030498241.2_17670__3483'],
            t_nested['Casgl158S09588__3522'],
            t_nested['CiLak.03G252300.1__32201'],
            t_nested['CiLak.04G185100.1__32201'],
            t_nested['lcl_CM025853.1_cds_KAB1213739.1_15127__262757'],
            t_nested['Qurub.09G170000.1__3512'],
            t_nested['lcl_NC_037093.1_cds_XP_024166822.1_38598__74649'],
            t_nested['Lj2g0027597.1__34305'],
            t_nested['Ca_17433__3827'],
            t_nested['Ler.1DRT.2g054060.1__41257'],
            t_nested['Medtr6g015950.1__3880'],
            t_nested['Vradi0190s00020.1__157791'],
            t_nested['Vradi04g09500.1__157791'],
            t_nested['Lj1g0018607.1__34305'],
            t_nested['Ca_20654__3827'],
            t_nested['lcl_OX451735.1_cds_CAI8591959.1_6550__3906'],
            t_nested['Vradi04g09490.1__157791'],
            t_nested['Ca_20655__3827'],
            t_nested['Medtr7g016090.2__3880'],
            t_nested['lcl_OX451735.1_cds_CAI8591957.1_6548__3906'],
            t_nested['Kaladp0046s0122.1.v1.1__63787'],
            t_nested['Lus10000711_PACid-23144352__4006'],
            t_nested['Soltu.DM.09G021170.1__4113'],
            t_nested['Lj1A65T52.1__105884'],
            t_nested['lcl_CM042060.1_cds_KAI3678813.1_38116__4217'],
            t_nested['AAM48133__F3H__flavonoid_pathway__137893'],
            t_nested['Lj1A65T51.1__105884'],
            t_nested['Pg_S2848.15__4054'],
            t_nested['lcl_NC_080145.1_cds_XP_057507794.1_23394__165200'],
            t_nested['Hma1.2p1_0273F.1_g109810.1__23110'],
            t_nested['pveT_jg6524.t1__170927'],
            t_nested['rna-gnl_WGS-JAMLDZ_EVM0012744.1__4058'],
            t_nested['Soltu.DM.06G028900.1__4113'],
            t_nested['rna-gnl_WGS-JAMLDZ_EVM0022992.1__4058'],
            t_nested['VIT_207s0005g03150.1__29760'],
            t_nested['lcl_NW_026775575.1_cds_XP_059663747.1_17704__4283'],
            t_nested['DCAR_005695__4039'],
            t_nested['HanXRQChr10g0281021__4232'],
            t_nested['lcl_JAMZMK010010835.1_cds_KAI7730131.1_28277__4212'],
            t_nested['rna-gnl_WGS-JAMLDZ_EVM0015450.1__4058'],
            t_nested['lcl_NC_080147.1_cds_XP_057512616.1_26605__165200'],
            t_nested['lcl_NW_026775575.1_cds_XP_059660363.1_17696__4283'],
            t_nested['Soltu.DM.02G013120.1__4113'],
            t_nested['lcl_CM061489.1_cds_KAK1426781.1_15461__13708'],
            t_nested['Pg_S0888.23__4054'],
            t_nested['DCAR_012540__4039'],
            t_nested['Lj2A1153T48.1__105884'],
            t_nested['Acc05348.1__3625'],
            t_nested['lcl_NC_039898.1_cds_XP_027114618.1_2067__13443'],
            t_nested['Oeu039683.1__4146'],
            t_nested['lcl_NC_032000.1_cds_XP_019256980.1_17445__49451'],
            t_nested['lcl_NC_062965.1_cds_XP_047977737.1_1129__49212'],
            t_nested['lcl_NW_026775575.1_cds_XP_059664398.1_17705__4283'],
            t_nested['pveT_jg25544.t1__170927'],
            t_nested['lcl_CM028329.1.g4555.t1__4392'],
            t_nested['Lj2A1153T49.1__105884'],
            t_nested['lcl_PKPP01011083.1_cds_PWA44819.1_54922__35608'],
            t_nested['DCAR_005694__4039'],
            t_nested['Pg_S2848.14__4054'],
            t_nested['DCAR_000938__4039'],
            t_nested['lcl_OX459122.1_cds_CAI9107661.1_20585__43536'],
            t_nested['lcl_NC_083345.1_cds_XP_060180555.1_37131__112863'],
            t_nested['lcl_OU503043.1_cds_CAI9766859.1_14290__56036'],
            t_nested['MSTRG.18317.1__82927'],
            t_nested['lcl_NC_062965.1_cds_XP_047961505.1_7237__49212'],
            t_nested['lcl_NW_026137595.1_cds_XP_051142067.1_22958__175694'],
            t_nested['lcl_NW_026137608.1_cds_XP_051150305.1_30360__175694'],
        ],
        [
            t_nested['Sspon.04G0019940-2D-mRNA-1__62335'],
            t_nested['BAA03647.1__IDS__mugineic_acid_biosynthesis__4513'],
        ],
        [
            t_nested['p5.00_sc00134_p0038.1__51953'],
            t_nested['CM031531.1.CM031531.1.g56644.t1__4682'],
            t_nested['evm.model.AsparagusV1_07.1877__4686'],
            t_nested['ASH_rna19788__1088818'],
        ],
        [
            t_nested['Vradi08g06870.1__157791'], # manually added nested 2ODD12B within 2ODD12
            t_nested['Medtr1g070135.1__3880'],
            t_nested['Ler.1DRT.1g071540.1__41257'],
            t_nested['Medtr1g070120.1__3880'],
            t_nested['Ler.1DRT.1g071570.1__41257'],
            t_nested['Medtr1g070080.1__3880'],
            t_nested['Ler.1DRT.1g071590.1__41257'],
            t_nested['Medtr1g070085.3__3880'],
        ],
        [
            t_nested['Cc01_g06970.1__49390'],
            t_nested['Cc01_g09870.1__49390'],
            t_nested['rna-gnl_WGS-JAMLDZ_EVM0004127.1__4058'],
        ],
        [
            t_nested['EF187826__H6H__scopolamine_biosynthesis__402998'],
            t_nested['AY356396.1__H6H__scopolamine_biosynthesis__243964'],
            t_nested['EF442802__H6H__scopolamine_biosynthesis__448155'],
            t_nested['EU530633.1__H6H__scopolamine_biosynthesis__512269'],
            t_nested['M62719__H6H__scopolamine_biosynthesis__4079'],
            t_nested['AF435417__H6H__scopolamine_biosynthesis__35625'],
        ]
    ]
    result = get_clusters(t_nested, parent_map=parent_map)
    assert result == expected_clusters_nestd_without_candidates

#%%

def test_get_clusters_nested_with_candidates():
    manual_nested_seqs_12B = [
        t_nested['Vradi08g06870.1__157791'],
        t_nested['Medtr1g070135.1__3880'],
        t_nested['Ler.1DRT.1g071540.1__41257'],
        t_nested['Medtr1g070120.1__3880'],
        t_nested['Ler.1DRT.1g071570.1__41257'],
        t_nested['Medtr1g070080.1__3880'],
        t_nested['Ler.1DRT.1g071590.1__41257'],
        t_nested['Medtr1g070085.3__3880'],
    ]

    manual_nested_seqs_12A = [
        t_nested['Cc01_g06970.1__49390'],
        t_nested['Cc01_g09870.1__49390'],
        t_nested['rna-gnl_WGS-JAMLDZ_EVM0004127.1__4058'],
    ]

    for leaf in t_nested:
        if leaf in manual_nested_seqs_12B:
            leaf.props["two_odd_id"] = "2ODD12B"
        elif leaf in manual_nested_seqs_12A:
            leaf.props["two_odd_id"] = "2ODD12A"

    manual_candidates = {
        'Sspon.04G0019940-2D-mRNA-1__62335': "2ODD12B",
        'p5.00_sc00134_p0038.1__51953': "2ODD12B",
        'ASH_rna19788__1088818': "2ODD12B",
        "Ler.1DRT.1g071570.1__41257": "2ODD12B",
        "Medtr1g070080.1__3880": "2ODD12B",
        "Soltu.DM.06G028900.1__4113": "2ODD12",
        'rna-gnl_WGS-JAMLDZ_EVM0022992.1__4058': "2ODD12",
        "EF187826__H6H__scopolamine_biosynthesis__402998": "2ODD12A",
        'EF442802__H6H__scopolamine_biosynthesis__448155': "2ODD12A",
        'lcl_CM059870.1_cds_KAK0581506.1_15868__4024': "2ODD12",
        'lcl_NW_006262339.1_cds_XP_024044614.1_12214__85681': "2ODD12",
        'Thecc.01G260000.1__3641': "2ODD12",
        'Dglom04181.1__34297': "2ODD12",
        'Blora11689.1__200023': "2ODD12",
        'Qurub.09G177400.1__3512': "2ODD12",
        'lcl_NC_083387.1_cds_XP_015891362.3_28814__326968': "2ODD12",
        'lcl_NC_083604.1_cds_XP_030498241.2_17670__3483': "2ODD12",
        'Casgl158S09588__3522': "2ODD12",
    }

    for leaf in t_nested:
        if leaf.name in manual_candidates:
            leaf.props["two_odd_id"] = "candidate"
    
    expected_clusters_nested_with_candidates = [

        [
            t_nested['Kaladp0089s0096.1.v1.1__63787'],
            t_nested['lcl_CM027381.1_cds_KAF9605627.1_22189__261450'],
            t_nested['Vadar_g6879.t1__229202'],
            t_nested['Acc25611.1__3625'],
            t_nested['Dglom04183.1__34297'],
            t_nested['Anaoc.0002s1042.1__171929'],
            t_nested['Eucgr.B01535.1__71139'],
            t_nested['Qurub.09G168900.1__3512'],
            t_nested['Prupe.8G225800.1__3760'],
            t_nested['Casgl267S13418__3522'],
            t_nested['Casgl528S19423__3522'],
            t_nested['lcl_CM025853.1_cds_KAB1213729.1_15117__262757'],
            t_nested['lcl_NC_083604.1_cds_XP_030496368.2_17710__3483'],
            t_nested['Lj5g0001478.1__34305'],
            t_nested['Lus10035782_PACid-23147456__4006'],
            t_nested['lcl_NW_011499926.1_cds_XP_011038934.1_22618__75702'],
            t_nested['Gobar.D08G228100.1__3634'],
            t_nested['Anaoc.1020s0004.1__171929'],
            t_nested['Ptrif.0003s3349.1__37690'],
            t_nested['lcl_CM059870.1_cds_KAK0580988.1_15865__4024'],
            t_nested['FvH4_2g21630.t4__57918'],
            t_nested['Sivu_ALN-8597__42043'],
            t_nested['Sila-35377__37657'],
            t_nested['Sc24g0003120.01__3999'],
            t_nested['TRINITY_DN229963_c3_g1_i1__122310'],
            t_nested['TRINITY_DN105553_c0_g1_i1__223224'],
            t_nested['Sc17g0000710.01__3999'],
            t_nested['DIACA2-20392__3570'],
            t_nested['Sivu_ALN-38268__42043'],
            t_nested['CM029399.1.CM029399.1.g1105.t1__3670'],
            t_nested['lcl_CM059870.1_cds_KAK0581506.1_15868__4024'],
            t_nested['lcl_NW_006262339.1_cds_XP_024044614.1_12214__85681'],
            t_nested['Thecc.01G260000.1__3641'],
            t_nested['Dglom04181.1__34297'],
            t_nested['Blora11689.1__200023'],
            t_nested['Qurub.09G177400.1__3512'],
            t_nested['lcl_NC_083387.1_cds_XP_015891362.3_28814__326968'],
            t_nested['lcl_NC_083604.1_cds_XP_030498241.2_17670__3483'],
            t_nested['Casgl158S09588__3522'],
            t_nested['CiLak.03G252300.1__32201'],
            t_nested['CiLak.04G185100.1__32201'],
            t_nested['lcl_CM025853.1_cds_KAB1213739.1_15127__262757'],
            t_nested['Qurub.09G170000.1__3512'],
            t_nested['lcl_NC_037093.1_cds_XP_024166822.1_38598__74649'],
            t_nested['Lj2g0027597.1__34305'],
            t_nested['Ca_17433__3827'],
            t_nested['Ler.1DRT.2g054060.1__41257'],
            t_nested['Medtr6g015950.1__3880'],
            t_nested['Vradi0190s00020.1__157791'],
            t_nested['Vradi04g09500.1__157791'],
            t_nested['Lj1g0018607.1__34305'],
            t_nested['Ca_20654__3827'],
            t_nested['lcl_OX451735.1_cds_CAI8591959.1_6550__3906'],
            t_nested['Vradi04g09490.1__157791'],
            t_nested['Ca_20655__3827'],
            t_nested['Medtr7g016090.2__3880'],
            t_nested['lcl_OX451735.1_cds_CAI8591957.1_6548__3906'],
            t_nested['Kaladp0046s0122.1.v1.1__63787'],
            t_nested['Lus10000711_PACid-23144352__4006'],
            t_nested['Soltu.DM.09G021170.1__4113'],
            t_nested['Lj1A65T52.1__105884'],
            t_nested['lcl_CM042060.1_cds_KAI3678813.1_38116__4217'],
            t_nested['AAM48133__F3H__flavonoid_pathway__137893'],
            t_nested['Lj1A65T51.1__105884'],
            t_nested['Pg_S2848.15__4054'],
            t_nested['lcl_NC_080145.1_cds_XP_057507794.1_23394__165200'],
            t_nested['Hma1.2p1_0273F.1_g109810.1__23110'],
            t_nested['pveT_jg6524.t1__170927'],
            t_nested['rna-gnl_WGS-JAMLDZ_EVM0012744.1__4058'],
            t_nested['rna-gnl_WGS-JAMLDZ_EVM0022992.1__4058'],
            t_nested['VIT_207s0005g03150.1__29760'],
            t_nested['lcl_NW_026775575.1_cds_XP_059663747.1_17704__4283'],
            t_nested['DCAR_005695__4039'],
            t_nested['HanXRQChr10g0281021__4232'],
            t_nested['lcl_JAMZMK010010835.1_cds_KAI7730131.1_28277__4212'],
            t_nested['rna-gnl_WGS-JAMLDZ_EVM0015450.1__4058'],
            t_nested['lcl_NC_080147.1_cds_XP_057512616.1_26605__165200'],
            t_nested['lcl_NW_026775575.1_cds_XP_059660363.1_17696__4283'],
            t_nested['Soltu.DM.02G013120.1__4113'],
            t_nested['lcl_CM061489.1_cds_KAK1426781.1_15461__13708'],
            t_nested['Pg_S0888.23__4054'],
            t_nested['DCAR_012540__4039'],
            t_nested['Lj2A1153T48.1__105884'],
            t_nested['Acc05348.1__3625'],
            t_nested['lcl_NC_039898.1_cds_XP_027114618.1_2067__13443'],
            t_nested['Oeu039683.1__4146'],
            t_nested['lcl_NC_032000.1_cds_XP_019256980.1_17445__49451'],
            t_nested['lcl_NC_062965.1_cds_XP_047977737.1_1129__49212'],
            t_nested['lcl_NW_026775575.1_cds_XP_059664398.1_17705__4283'],
            t_nested['pveT_jg25544.t1__170927'],
            t_nested['lcl_CM028329.1.g4555.t1__4392'],
            t_nested['Lj2A1153T49.1__105884'],
            t_nested['lcl_PKPP01011083.1_cds_PWA44819.1_54922__35608'],
            t_nested['DCAR_005694__4039'],
            t_nested['Pg_S2848.14__4054'],
            t_nested['DCAR_000938__4039'],
            t_nested['lcl_OX459122.1_cds_CAI9107661.1_20585__43536'],
            t_nested['lcl_NC_083345.1_cds_XP_060180555.1_37131__112863'],
            t_nested['lcl_OU503043.1_cds_CAI9766859.1_14290__56036'],
            t_nested['MSTRG.18317.1__82927'],
            t_nested['lcl_NC_062965.1_cds_XP_047961505.1_7237__49212'],
            t_nested['lcl_NW_026137595.1_cds_XP_051142067.1_22958__175694'],
            t_nested['lcl_NW_026137608.1_cds_XP_051150305.1_30360__175694'],
        ],
        [
            t_nested['Sspon.04G0019940-2D-mRNA-1__62335'],
            t_nested['BAA03647.1__IDS__mugineic_acid_biosynthesis__4513'],
        ],
        [
            t_nested['p5.00_sc00134_p0038.1__51953'],
            t_nested['CM031531.1.CM031531.1.g56644.t1__4682'],
            t_nested['evm.model.AsparagusV1_07.1877__4686'],
            t_nested['ASH_rna19788__1088818'],
        ],
        [
            t_nested['Vradi08g06870.1__157791'], # manually added nested 2ODD12B within 2ODD12
            t_nested['Medtr1g070135.1__3880'],
            t_nested['Ler.1DRT.1g071540.1__41257'],
            t_nested['Medtr1g070120.1__3880'],
            t_nested['Ler.1DRT.1g071570.1__41257'],
            t_nested['Medtr1g070080.1__3880'],
            t_nested['Ler.1DRT.1g071590.1__41257'],
            t_nested['Medtr1g070085.3__3880'],
        ],
        [
            t_nested['Cc01_g06970.1__49390'],
            t_nested['Cc01_g09870.1__49390'],
            t_nested['rna-gnl_WGS-JAMLDZ_EVM0004127.1__4058'],
        ],
        [
            t_nested['Soltu.DM.06G028900.1__4113'],
            t_nested['EF187826__H6H__scopolamine_biosynthesis__402998'],
            t_nested['AY356396.1__H6H__scopolamine_biosynthesis__243964'],
            t_nested['EF442802__H6H__scopolamine_biosynthesis__448155'],
            t_nested['EU530633.1__H6H__scopolamine_biosynthesis__512269'],
            t_nested['M62719__H6H__scopolamine_biosynthesis__4079'],
            t_nested['AF435417__H6H__scopolamine_biosynthesis__35625'],
        ]
    ]

    result = get_clusters(t_nested, dist_dict=dist_dict_nested, parent_map=parent_map)
    assert result == expected_clusters_nested_with_candidates


# %% ========================== Test compute_cluster_neighbors ==========================
dist_dict = build_distance_lookup(t)

#get median distance between all pairs of leaves in the tree
all_distances = []
for leaf1 in t.leaves():
    for leaf2 in t.leaves():
        if leaf1 != leaf2:
            all_distances.append(dist_dict[leaf1.name][leaf2.name])
median_distance = sorted(all_distances)[len(all_distances) // 2]



def test_compute_cluster_neighbors():

    expected = {'0': {'closest_cluster': 1, 'distance': 1.708},
                '1': {'closest_cluster': 28, 'distance': 1.672},
                '2': {'closest_cluster': 5, 'distance': 1.572},
                '3': {'closest_cluster': 5, 'distance': 1.656},
                '4': {'closest_cluster': 5, 'distance': 1.318},
                '5': {'closest_cluster': 6, 'distance': 1.131},
                '6': {'closest_cluster': 5, 'distance': 1.131},
                '7': {'closest_cluster': 6, 'distance': 1.493},
                '8': {'closest_cluster': 9, 'distance': 1.443},
                '9': {'closest_cluster': 14, 'distance': 1.015},
                '10': {'closest_cluster': 14, 'distance': 1.105},
                '11': {'closest_cluster': 14, 'distance': 0.817},
                '12': {'closest_cluster': 14, 'distance': 0.79},
                '13': {'closest_cluster': 14, 'distance': 0.951},
                '14': {'closest_cluster': 12, 'distance': 0.79},
                '15': {'closest_cluster': 14, 'distance': 0.81},
                '16': {'closest_cluster': 14, 'distance': 1.075},
                '17': {'closest_cluster': 14, 'distance': 1.148},
                '18': {'closest_cluster': 14, 'distance': 0.938},
                '19': {'closest_cluster': 14, 'distance': 0.938},
                '20': {'closest_cluster': 14, 'distance': 1.073},
                '21': {'closest_cluster': 14, 'distance': 1.176},
                '22': {'closest_cluster': 14, 'distance': 1.0},
                '23': {'closest_cluster': 14, 'distance': 1.17},
                '24': {'closest_cluster': 25, 'distance': 1.132},
                '25': {'closest_cluster': 14, 'distance': 1.009},
                '26': {'closest_cluster': 27, 'distance': 0.987},
                '27': {'closest_cluster': 26, 'distance': 0.987},
                '28': {'closest_cluster': 9, 'distance': 1.177}}
    
    manual_candidates = {
    "lcl_NC_084852.1_cds_XP_061960601.1_3168__3691" : "2ODD34",
    'VIT_204s0008g04920.2__29760' : '2ODD34',
    'Ptrif.0001s0318.1__37690' : '2ODD34', 
    'Atru_chr7_2342__47965' : '2ODD34',
    'lcl_NW_019168159.1_cds_XP_022728803.1_47147__66656' : '2ODD34',
    'lcl_NC_045134.1_cds_XP_031407529.1_33286__22663' : '2ODD34',
    'GWHTACBH000860__413952' : '2ODD34',
    'Potri.006G248000.1__3694' : '2ODD34',
    'lcl_NC_065570.1_cds_XP_050209042.1_2787__3986' : '2ODD34',
    'Lus10008097_PACid-23169790__4006' : '2ODD34', 
    'Acc09929.1__3625': "minor_2ODD_cluster",
    'lcl_NC_031989.1_cds_XP_019252277.1_2282__49451': "2ODD35",
    'Kaladp0011s0492.1.v1.1__63787': "2ODD35",
    'lcl_OX459121.1_cds_CAI9101307.1_14231__43536': "2ODD35",
    'rna-gnl_WGS-JAMLDZ_EVM0034906.1__4058': "2ODD35",
    'lcl_NC_039901.1_cds_XP_027112943.1_15453__13443': "2ODD35",
    'rna-gnl_WGS-JAMLDZ_EVM0007978.1__4058': "2ODD35", 
    'Bol038153__3712': "2ODD36",
    "lcl_NC_053024.1_cds_XP_048567703.1_14633__4572": "unresolved"
    }

    for leaf in t:
        if leaf.name in manual_candidates:
            leaf.props["two_odd_id"] = "candidate"

    clusters = get_clusters(t)

    result = compute_cluster_neighbors(clusters=clusters,
                                       dist_dict=dist_dict)
    assert result == expected

#%% 


def test_two_odd_id_to_landscape_indices():
    manual_candidates = {
        "lcl_NC_084852.1_cds_XP_061960601.1_3168__3691" : "2ODD34",
        'VIT_204s0008g04920.2__29760' : '2ODD34',
        'Ptrif.0001s0318.1__37690' : '2ODD34',
        'Atru_chr7_2342__47965' : '2ODD34',
        'lcl_NW_019168159.1_cds_XP_022728803.1_47147__66656' : '2ODD34',
        'lcl_NC_045134.1_cds_XP_031407529.1_33286__22663' : '2ODD34',
        'GWHTACBH000860__413952' : '2ODD34',
        'Potri.006G248000.1__3694' : '2ODD34',
        'lcl_NC_065570.1_cds_XP_050209042.1_2787__3986' : '2ODD34',
        'Lus10008097_PACid-23169790__4006' : '2ODD34', 
        'Acc09929.1__3625': "minor_2ODD_cluster",
        'lcl_NC_031989.1_cds_XP_019252277.1_2282__49451': "2ODD35",
        'Kaladp0011s0492.1.v1.1__63787': "2ODD35",
        'lcl_OX459121.1_cds_CAI9101307.1_14231__43536': "2ODD35",
        'rna-gnl_WGS-JAMLDZ_EVM0034906.1__4058': "2ODD35",
        'lcl_NC_039901.1_cds_XP_027112943.1_15453__13443': "2ODD35",
        'rna-gnl_WGS-JAMLDZ_EVM0007978.1__4058': "2ODD35", 
        'Bol038153__3712': "2ODD36",
        "lcl_NC_053024.1_cds_XP_048567703.1_14633__4572": "unresolved"
        }

    for leaf in t:
        if leaf.name in manual_candidates:
            leaf.props["two_odd_id"] = "candidate"

    clusters = get_clusters(t)

    expected = {'candidates_only': [0],
                '2ODD34': [1, 2, 3, 4, 5, 6],
                'minor_2ODD_cluster': [7, 11, 15, 17, 18, 21, 22, 23],
                '2ODD38': [8, 9, 10, 28],
                '2ODD37': [12, 13, 14, 16, 19],
                '2ODD35': [20, 27],
                '2ODD36': [24, 25, 26]}

    result = two_odd_id_to_cluster_indices(clusters)
    assert result == expected


def test_seq_to_cluster_idx():

    manual_candidates = {
        "lcl_NC_084852.1_cds_XP_061960601.1_3168__3691" : "2ODD34",
        'VIT_204s0008g04920.2__29760' : '2ODD34',
        'Ptrif.0001s0318.1__37690' : '2ODD34',
        'Atru_chr7_2342__47965' : '2ODD34',
        'lcl_NW_019168159.1_cds_XP_022728803.1_47147__66656' : '2ODD34',
        'lcl_NC_045134.1_cds_XP_031407529.1_33286__22663' : '2ODD34',
        'GWHTACBH000860__413952' : '2ODD34',
        'Potri.006G248000.1__3694' : '2ODD34',
        'lcl_NC_065570.1_cds_XP_050209042.1_2787__3986' : '2ODD34',
        'Lus10008097_PACid-23169790__4006' : '2ODD34', 
        'Acc09929.1__3625': "minor_2ODD_cluster",
        'lcl_NC_031989.1_cds_XP_019252277.1_2282__49451': "2ODD35",
        'Kaladp0011s0492.1.v1.1__63787': "2ODD35",
        'lcl_OX459121.1_cds_CAI9101307.1_14231__43536': "2ODD35",
        'rna-gnl_WGS-JAMLDZ_EVM0034906.1__4058': "2ODD35",
        'lcl_NC_039901.1_cds_XP_027112943.1_15453__13443': "2ODD35",
        'rna-gnl_WGS-JAMLDZ_EVM0007978.1__4058': "2ODD35", 
        'Bol038153__3712': "2ODD36",
        "lcl_NC_053024.1_cds_XP_048567703.1_14633__4572": "unresolved"
        }

    for leaf in t:
        if leaf.name in manual_candidates:
            leaf.props["two_odd_id"] = "candidate"

    clusters = get_clusters(t)

    expected = {'lcl_NC_084852.1_cds_XP_061960601.1_3168__3691': 0,
                'rna-XM_031626300.2__210225': 1,
                'CKAN_00595100__337451': 2,
                'Aristolochia_fimbriatalcl_CM034078.1_cds_KAG9459387.1_3895__158543': 2,
                'Gobar.D13G239800.1__3634': 2,
                'Lj8A611G43.1__105884': 2,
                'FvH4_2g29990.t1__57918': 2,
                'Casgl23S04257__3522': 2,
                'CiLak.07G220400.1__32201': 2,
                'lcl_NC_065570.1_cds_XP_050209022.1_2785__3986': 2,
                'Lus10013130_PACid-23164132__4006': 2,
                'evm_27.model.AmTr_v1.0_scaffold00174.4__13333': 3,
                'RZC71586__3469': 4,
                'lcl_NW_010729080.1_cds_XP_010249917.1_5454__4432': 5,
                'rna-XM_042628350.1__60698': 5,
                'lcl_CM056811.1_cds_KAJ8637692.1_11945__3435': 5,
                'Aristolochia_fimbriatalcl_CM034078.1_cds_KAG9459386.1_3894__158543': 5,
                'Kaladp1006s0012.1.v1.1__63787': 6,
                'VIT_204s0008g04920.2__29760': 6,
                'Ptrif.0001s0318.1__37690': 6,
                'Atru_chr7_2342__47965': 6,
                'lcl_NW_019168159.1_cds_XP_022728803.1_47147__66656': 6,
                'lcl_NC_045134.1_cds_XP_031407529.1_33286__22663': 6,
                'GWHTACBH000860__413952': 6,
                'Potri.006G248000.1__3694': 6,
                'lcl_NC_065570.1_cds_XP_050209042.1_2787__3986': 6,
                'Lus10008097_PACid-23169790__4006': 6,
                'Acc06446.1__3625': 6,
                'Acc16585.1__3625': 6,
                'pveT_jg24143.t1__170927': 6,
                'pveT_jg28701.t1__170927': 6,
                'lcl_NC_083383.1_cds_XP_048332159.2_17459__326968': 6,
                'Cnepa31456.1__79760': 6,
                'Umino40998.1__262084': 6,
                'FvH4_2g29991.t1__57918': 6,
                'Lj1g0022518.1__34305': 6,
                'Lj4g0022664.1__34305': 6,
                'Medtr1g011600.1__3880': 6,
                'Ca_10933__3827': 6,
                'Ler.1DRT.1g009880.1__41257': 6,
                'Ler.1DRT.3g073930.1__41257': 6,
                'Medtr3g108520.1__3880': 6,
                'Ca_23071__3827': 6,
                'Ca_01791__3827': 6,
                'Medtr4g021360.1__3880': 6,
                'Ca_05192__3827': 6,
                'Medtr4g021380.1__3880': 6,
                'lcl_OX451737.1_cds_CAI8600332.1_14923__3906': 6,
                'pveT_jg20218.t1__170927': 7,
                'pveT_jg20220.t1__170927': 7,
                'Acc09929.1__3625': 7,
                'pveT_jg33784.t1__170927': 7,
                'Eucgr.C03724.1__71139': 8,
                'lcl_NC_045130.1_cds_XP_031390300.1_19564__22663': 8,
                'Eucgr.C04147.1__71139': 8,
                'Acc23206.1__3625': 9,
                'lcl_CM027379.1_cds_KAF9612985.1_11655__261450': 9,
                'RZC71605__3469': 9,
                'rna-XM_042647017.1__60698': 9,
                'lcl_CM027379.1_cds_KAF9612671.1_9820__261450': 10,
                'Lj7A541T77.1__105884': 11,
                'Pg_S7309.1__4054': 12,
                'DCAR_015859__4039': 12,
                'CsatW809539.1__3659': 13,
                'lcl_NW_026775580.1_cds_XP_059633894.1_39517__4283': 14,
                'Dglom27812.1__34297': 14,
                'lcl_CM059866.1_cds_KAK0595785.1_1317__4024': 14,
                'Ptrif.0007s1329.1__37690': 14,
                'Potri.001G176500.1__3694': 14,
                'lcl_NC_065574.1_cds_XP_050234707.1_20560__3986': 14,
                'Umino20162.1__262084': 14,
                'lcl_NC_083379.1_cds_XP_015868728.2_4750__326968': 14,
                'lcl_NC_084800.1_cds_XP_062086870.1_35676__3486': 14,
                'lcl_NC_045131.1_cds_XP_031396687.1_21997__22663': 14,
                'Eucgr.J03022.1__71139': 14,
                'VIT_209s0002g05290.1__29760': 14,
                'lcl_NW_026775580.1_cds_XP_059636056.1_39514__4283': 14,
                'Acc08999.1__3625': 14,
                'Acc07679.1__3625': 14,
                'FSB010835101__28930': 14,
                'DCAR_015858__4039': 14,
                'Casgl236S24615__3522': 15,
                'CiLak.15G093300.1__32201': 15,
                'Qurub.10G093800.1__3512': 15,
                'lcl_CM025849.1_cds_KAB1226463.1_1752__262757': 15,
                'Tcord05348.2__703396': 15,
                'HU01G00826.1__176265': 15,
                'Sc13g0006470.01__3999': 15,
                'TRINITY_DN108901_c0_g1_i1__223224': 15,
                'CM008280.1.CM008280.1.g987.t1__62330': 15,
                'pveT_jg9723.t1__170927': 16,
                'Cc06_g15940.1__49390': 17,
                'Gobar.D05G137100.1__3634': 18,
                'lcl_OX459118.1_cds_CAI9089985.1_2909__43536': 18,
                'TRINITY_DN197356_c3_g1_i1__4102': 19,
                'lcl_CM061488.1_cds_KAK1430368.1_13046__13708': 19,
                'Oeu000043.1__4146': 19,
                'lcl_NW_026137546.1_cds_XP_051123007.1_5823__175694': 19,
                'lcl_NW_025952766.1_cds_XP_047954044.1_46391__49212': 19,
                'lcl_NC_028637.1_cds_XP_015064474.1_2338__28526': 19,
                'Cc00_g25210.1__49390': 19,
                'rna-gnl_WGS-JAMLDZ_EVM0036354.1__4058': 19,
                'lcl_OX459121.1_cds_CAI9101309.1_14233__43536': 19,
                'lcl_NC_031989.1_cds_XP_019252277.1_2282__49451': 20,
                'Kaladp0011s0492.1.v1.1__63787': 20,
                'lcl_OX459121.1_cds_CAI9101307.1_14231__43536': 20,
                'rna-gnl_WGS-JAMLDZ_EVM0034906.1__4058': 20,
                'lcl_NC_039901.1_cds_XP_027112943.1_15453__13443': 20,
                'rna-gnl_WGS-JAMLDZ_EVM0007978.1__4058': 20,
                'lcl_NC_080155.1_cds_XP_057473373.1_43266__165200': 20,
                'Hma1.2p1_1763F.1_g301360.1__23110': 20,
                'pveT_jg21836.t1__170927': 20,
                'Lj5A64G50.1__105884': 20,
                'lcl_CM028326.1.g5350.t1__4392': 20,
                'rna-gnl_WGS-JAMLDZ_EVM0017794.1__4058': 20,
                'lcl_NC_039905.1_cds_XP_027126787.1_26673__13443': 20,
                'lcl_OX459124.1_cds_CAI9114063.1_26987__43536': 20,
                'lcl_PKPP01007341.1_cds_PWA53652.1_46097__35608': 20,
                'rna-gnl_WGS-JAMLDZ_EVM0014154.1__4058': 20,
                'lcl_OX459124.1_cds_CAI9114065.1_26989__43536': 20,
                'Soltu.DM.12G010750.2__4113': 20,
                'lcl_OU503036.1_cds_CAI9753676.1_1107__56036': 20,
                'lcl_NC_062966.1_cds_XP_047959520.1_10496__49212': 20,
                'lcl_NW_026137468.1_cds_XP_051115997.1_3322__175694': 20,
                'Cc02_g31130.1__49390': 21,
                'Cc00_g28750.1__49390': 21,
                'rna-gnl_WGS-JAMLDZ_EVM0024530.1__4058': 21,
                'Cc03_g08070.1__49390': 21,
                'Anaoc.0007s0766.1__171929': 22,
                'Atru_chr1_3170__47965': 22,
                'FSB015489301__28930': 22,
                'lcl_CM059866.1_cds_KAK0596815.1_1321__4024': 22,
                'FSB011836701__28930': 22,
                'CiLak.01G168500.1__32201': 22,
                'Qurub.10G026200.1__3512': 22,
                'Qurub.10G071600.1__3512': 22,
                'lcl_CM025853.1_cds_KAB1214321.1_15709__262757': 22,
                'CiLak.01G168600.1__32201': 22,
                'lcl_NC_044913.1_cds_XP_030940120.1_44897__97700': 22,
                'Prupe.5G048100.1__3760': 22,
                'Umino00085.1__262084': 22,
                'lcl_NC_024133.1_cds_XP_008243039.1_26638__102107': 22,
                'lcl_NC_084796.1_cds_XP_062119410.1_20878__3486': 23,
                'Umino11360.1__262084': 23,
                'FvH4_3g36530.t1__57918': 23,
                'lcl_NC_083380.1_cds_XP_015881543.2_5733__326968': 23,
                'Kaladp0007s0029.1.v1.1__63787': 23,
                'Blora12794.1__200023': 23,
                'VIT_209s0002g05340.1__29760': 24,
                'lcl_NW_017353139.1_cds_XP_018486826.1_3787__3726': 24,
                'lcl_NC_045129.1_cds_XP_031385398.1_13438__22663': 24,
                'lcl_JAGHRR010000128.1_cds_KAI3418995.1_3930__120290': 24,
                'lcl_NW_026137602.1_cds_XP_051147854.1_28335__175694': 24,
                'Hma1.2p1_0515F.1_g170445.1__23110': 24,
                'lcl_CM028334.1.g444.t1__4392': 24,
                'lcl_OU503038.1_cds_CAI9757760.1_5191__56036': 24,
                'Vadar_g3125.t1__229202': 24,
                'Acc17837.1__3625': 24,
                'pveT_jg21833.t1__170927': 24,
                'Lj5A64T54.1__105884': 24,
                'lcl_NW_026775581.1_cds_XP_059644652.1_44213__4283': 24,
                'Cc03_g10170.1__49390': 24,
                'lcl_CM061492.1_cds_KAK1416820.1_25937__13708': 24,
                'DCAR_016236__4039': 24,
                'lcl_NC_062967.1_cds_XP_047968971.1_19558__49212': 24,
                'TRINITY_DN176594_c2_g1_i2__4102': 24,
                'lcl_NC_045131.1_cds_XP_031399179.1_24312__22663': 24,
                'Eucgr.F02562.1__71139': 24,
                'lcl_NW_006262075.1_cds_XP_006427129.2_25983__85681': 24,
                'lcl_CM059866.1_cds_KAK0598226.1_1319__4024': 24,
                'lcl_NC_052253.1_cds_XP_038693680.1_44458__458696': 24,
                'lcl_NW_011501060.1_cds_XP_011013043.1_5386__75702': 24,
                'lcl_NC_045129.1_cds_XP_031385404.1_13435__22663': 24,
                'Eucgr.F02566.1__71139': 24,
                'Bol038153__3712': 24,
                'Lj1g0016927.1__34305': 24,
                'Ler.1DRT.6g050030.1__41257': 24,
                'Medtr7g090520.1__3880': 24,
                'Vradi08g02800.1__157791': 24,
                'Lj1g0005382.1__34305': 24,
                'Vradi08g02780.1__157791': 24,
                'lcl_CM059866.1_cds_KAK0596259.1_1320__4024': 24,
                'Ptrif.0007s1335.1__37690': 24,
                'Anaoc.0007s0767.1__171929': 24,
                'Gorai.002G057700.1_PACid-26792653__29730': 24,
                'Ptrif.0007s1331.1__37690': 24,
                'lcl_NC_052243.1_cds_XP_038718244.1_22960__458696': 24,
                'lcl_NC_083379.1_cds_XP_015868729.2_4751__326968': 24,
                'lcl_NC_037092.1_cds_XP_024158425.1_33959__74649': 24,
                'Umino20160.1__262084': 24,
                'Umino20144.1__262084': 24,
                'lcl_NC_083601.1_cds_XP_030488467.2_4391__3483': 24,
                'Pparo25676.1__386216': 24,
                'Dglom27811.1__34297': 24,
                'Blora12130.2__200023': 24,
                'lcl_NC_044908.1_cds_XP_030969963.1_22137__97700': 24,
                'lcl_NC_052253.1_cds_XP_038693679.1_44455__458696': 24,
                'Lus10021025_PACid-23182311__4006': 24,
                'lcl_CM058159.1_cds_KAJ9692491.1_17293__103349': 25,
                'Vradi08g02810.1__157791': 25,
                'Ler.1DRT.1g027930.1__41257': 25,
                'Lj5g0008310.1__34305': 25,
                'TRINITY_DN114504_c0_g1_i1__4102': 25,
                'Eucgr.J03023.1__71139': 25,
                'lcl_NC_045131.1_cds_XP_031396688.1_21996__22663': 25,
                'Gorai.011G079100.1_PACid-26810489__29730': 25,
                'Thhalv10011669m__72664': 25,
                'Anaoc.0007s0768.1__171929': 25,
                'lcl_NW_006262075.1_cds_XP_006427132.2_25985__85681': 25,
                'Atru_chr1_3180__47965': 25,
                'Atru_chr1_3175__47965': 25,
                'Potri.001G176000.1__3694': 25,
                'lcl_NC_052244.1_cds_XP_038721828.1_25425__458696': 25,
                'lcl_NC_065574.1_cds_XP_050234356.1_20365__3986': 25,
                'lcl_NC_044908.1_cds_XP_030968259.1_22138__97700': 25,
                'lcl_CM025851.1_cds_KAB1218619.1_7986__262757': 25,
                'CiLak.01G242800.1__32201': 25,
                'lcl_NC_037089.1_cds_XP_024183820.1_11602__74649': 25,
                'lcl_NC_083379.1_cds_XP_048325586.2_4752__326968': 25,
                'lcl_NC_084800.1_cds_XP_062086877.1_35687__3486': 25,
                'Umino20121.1__262084': 25,
                'Dglom27810.1__34297': 25,
                'Tcord01333.1__703396': 25,
                'Blora12127.1__200023': 25,
                'CiLak.02G150500.1__32201': 26,
                'FSB014894201__28930': 26,
                'lcl_NC_044908.1_cds_XP_030968569.1_22136__97700': 26,
                'lcl_PKPP01002761.1_cds_PWA73280.1_26522__35608': 27,
                'rna-gnl_WGS-JAMLDZ_EVM0012381.1__4058': 27,
                'lcl_CM028334.1.g446.t1__4392': 27,
                'Lj5A64T53.1__105884': 27,
                'Hma1.2p1_3098F.1_g332680.1__23110': 27,
                'lcl_NW_026775581.1_cds_XP_059639306.1_44215__4283': 27,
                'Vadar_g3126.t1__229202': 27,
                'Acc12585.1__3625': 27,
                'pveT_jg21835.t1__170927': 27,
                'Pg_S0988.1__4054': 27,
                'DCAR_016235__4039': 27,
                'lcl_OX459124.1_cds_CAI9114062.1_26986__43536': 27,
                'lcl_OU503036.1_cds_CAI9753675.1_1106__56036': 27,
                'lcl_NC_062966.1_cds_XP_047959442.1_10495__49212': 27,
                'lcl_NW_026137602.1_cds_XP_051147853.1_28334__175694': 27,
                'DCAR_007548__4039': 27,
                'lcl_CM061487.1_cds_KAK1434086.1_11004__13708': 27,
                'HanXRQChr09g0275101__4232': 27,
                'Lj5C65T2.1__105884': 27,
                'lcl_NW_026775581.1_cds_XP_059644765.1_44212__4283': 27,
                'lcl_NC_080143.1_cds_XP_057503960.1_20127__165200': 27,
                'pveT_jg21832.t1__170927': 27,
                'lcl_CM028334.1.g443.t1__4392': 27,
                'rna-gnl_WGS-JAMLDZ_EVM0018858.1__4058': 27,
                'SMEL_006g260960.1.01__4111': 27,
                'lcl_NC_039905.1_cds_XP_027126923.1_26675__13443': 27,
                'lcl_OU503036.1_cds_CAI9753673.1_1104__56036': 27,
                'MSTRG.22708.1__82927': 27,
                'lcl_NC_062966.1_cds_XP_047960078.1_10492__49212': 27,
                'lcl_NW_026137587.1_cds_XP_051135949.1_16701__175694': 27,
                'lcl_NW_026137587.1_cds_XP_051135915.1_16702__175694': 27,
                'lcl_CM008454.1_cds_PHT32940.1_29277__33114': 28,
                'lcl_CM008444.1_cds_PHT55775.1_4261__33114': 28,
                'TRINITY_DN232929_c4_g6_i3__4102': 28,
                'lcl_CM008445.1_cds_PHT51893.1_6355__33114': 28,
                'lcl_NC_083337.1_cds_XP_060201071.1_373__112863': 28,
                'lcl_NC_083338.1_cds_XP_060200949.1_5353__112863': 28,
                'Soltu.DM.01G002190.1__4113': 28,
                'SMEL_001g152420.1.01__4111': 28,
                'lcl_NC_028638.1_cds_XP_015067085.1_5485__28526': 28,
                'TRINITY_DN725607_c0_g1_i1__189803': 28,
                'Vadar_g43141.t1__229202': 28,
                'KT390173__DPS__etoposide_biosynthesis__93608': 28,
                'lcl_NW_026775571.1_cds_XP_059664661.1_3040__4283': 28,
                'AF417859__AOP3__glucosinolate_biosynthesis__3702': 28,
                'AF417858__AOP2__glucosinolate_biosynthesis__3702': 28,
                'lcl_CM042048.1_cds_KAI3759661.1_7640__4217': 28,
                'lcl_NW_026137594.1_cds_XP_051141059.1_21911__175694': 28,
                'lcl_NC_083379.1_cds_XP_048320867.2_1719__326968': 28,
                'Umino20164.1__262084': 28,
                'Cc08_g04750.1__49390': 28,
                'Lus10023024_PACid-23146145__4006': 28,
                'lcl_NC_045131.1_cds_XP_031396066.1_23250__22663': 28,
                'lcl_NC_052241.1_cds_XP_038713292.1_19041__458696': 28,
                'Gobar.A11G255200.1__3634': 28,
                'Medtr4g011690.1__3880': 28,
                'Ca_13466__3827': 28,
                'Ler.1DRT.4g003660.1__41257': 28,
                'Ca_13467__3827': 28,
                'Lj3g0013768.1__34305': 28,
                'lcl_OX451740.1_cds_CAI8613451.1_28042__3906': 28,
                'Lj1g0004205.1__34305': 28,
                'Vradi0154s00060.1__157791': 28,
                'lcl_NC_052245.1_cds_XP_038722105.1_29001__458696': 28,
                'lcl_NC_065574.1_cds_XP_050231702.1_19581__3986': 28,
                'Lus10005516_PACid-23158832__4006': 28,
                'lcl_NC_052393.1_cds_XP_008784635.1_5049__42345': 28,
                'Vadar_g3127.t1__229202': 28,
                'VIT_209s0002g05280.1__29760': 28,
                'Acc17838.1__3625': 28,
                'lcl_NW_026775581.1_cds_XP_059643089.1_44218__4283': 28,
                'lcl_CM028334.1.g447.t1__4392': 28,
                'Lj9C505T3.1__105884': 28,
                'Pg_S0213.71__4054': 28,
                'DCAR_016234__4039': 28,
                'lcl_CM027379.1_cds_KAF9612984.1_11656__261450': 28,
                'lcl_PKPP01007341.1_cds_PWA53651.1_46096__35608': 28,
                'Hma1.2p1_0515F.1_g170455.1__23110': 28,
                'Ah1G51780.1__81970': 28,
                'Cvari18296.1__869952': 28,
                'Umino20166.1__262084': 28,
                'lcl_NC_083379.1_cds_XP_015868730.2_4749__326968': 28,
                'lcl_NC_044908.1_cds_XP_030973050.1_22133__97700': 28,
                'CiLak.01G242700.1__32201': 28,
                'lcl_CM025851.1_cds_KAB1218616.1_7983__262757': 28,
                'Prupe.3G075700.1__3760': 28,
                'Dglom27813.1__34297': 28,
                'Blora09835.1__200023': 28,
                'lcl_NW_006262075.1_cds_XP_006427127.2_25981__85681': 28,
                'Potri.001G176200.1__3694': 28,
                'lcl_NC_065574.1_cds_XP_050231462.1_20561__3986': 28,
                'lcl_JAGHRR010000039.1_cds_KAI3437248.1_15716__120290': 28,
                'Anaoc.0007s0764.1__171929': 28,
                'lcl_CM059866.1_cds_KAK0595174.1_1315__4024': 28,
                'OMO54201__210143': 28,
                'lcl_NC_052253.1_cds_XP_038694945.1_44457__458696': 28,
                'Lus10016659_PACid-23143917__4006': 28,
                'pveT_jg11526.t1__170927': 28,
                'lcl_NC_083601.1_cds_XP_060957775.1_4387__3483': 28,
                'lcl_OX459123.1_cds_CAI9108625.1_21549__43536': 28,
                'Kaladp0008s0237.1.v1.1__63787': 28,
                'Soltu.DM.06G023510.1__4113': 28,
                'rna-gnl_WGS-JAMLDZ_EVM0014967.1__4058': 28,
                'lcl_OU503036.1_cds_CAI9753677.1_1108__56036': 28,
                'MSTRG.11101.1__82927': 28,
                'lcl_NC_062966.1_cds_XP_047958340.1_10499__49212': 28,
                'lcl_NC_062965.1_cds_XP_047956864.1_861__49212': 28,
                'Mba05_g08120.1.v1.1__52838': 28,
                'CALSI_Maker00003989__746888': 28,
                'Mba03_g21520.1.v1.1__52838': 28,
                'evm.model.AsparagusV1_08.662__4686': 28,
                'Ma04_t25080.1__4641': 28,
                'lcl_CP136891.1_cds_WOK98631.1_7322__4628': 28,
                'Mba05_g08100.1.v1.1__52838': 28,
                'Mtr04_g32474.1__320322': 28,
                'lcl_CP136893.1_cds_WOL03332.1_12022__4628': 28,
                'CALSI_Maker00009131__746888': 28,
                'ORUFI05G28750.1__4529': 28,
                'ORUFI11G09650.1__4529': 28,
                'lcl_NC_052396.1_cds_XP_008794218.1_12468__42345': 28,
                'Spipo4G0051700__29656': 28,
                'lcl_CP136890.1_cds_WOK92059.1_750__4628': 28,
                'CALSI_Maker00050580__746888': 28,
                'evm.model.AsparagusV1_01.591__4686': 28,
                'ASH_rna4294__1088818': 28,
                'strangu_020525-RA__38733': 28,
                'lcl_NC_053024.1_cds_XP_048567703.1_14633__4572': 28}
    result = seq_to_cluster_idx(clusters)
    assert result == expected 



def test_cluster_meta_info():
    manual_candidates = {
        "lcl_NC_084852.1_cds_XP_061960601.1_3168__3691" : "2ODD34",
        'VIT_204s0008g04920.2__29760' : '2ODD34',
        'Ptrif.0001s0318.1__37690' : '2ODD34',
        'Atru_chr7_2342__47965' : '2ODD34',
        'lcl_NW_019168159.1_cds_XP_022728803.1_47147__66656' : '2ODD34',
        'lcl_NC_045134.1_cds_XP_031407529.1_33286__22663' : '2ODD34',
        'GWHTACBH000860__413952' : '2ODD34',
        'Potri.006G248000.1__3694' : '2ODD34',
        'lcl_NC_065570.1_cds_XP_050209042.1_2787__3986' : '2ODD34',
        'Lus10008097_PACid-23169790__4006' : '2ODD34', 
        'Acc09929.1__3625': "minor_2ODD_cluster",
        'lcl_NC_031989.1_cds_XP_019252277.1_2282__49451': "2ODD35",
        'Kaladp0011s0492.1.v1.1__63787': "2ODD35",
        'lcl_OX459121.1_cds_CAI9101307.1_14231__43536': "2ODD35",
        'rna-gnl_WGS-JAMLDZ_EVM0034906.1__4058': "2ODD35",
        'lcl_NC_039901.1_cds_XP_027112943.1_15453__13443': "2ODD35",
        'rna-gnl_WGS-JAMLDZ_EVM0007978.1__4058': "2ODD35", 
        'Bol038153__3712': "2ODD36",
        "lcl_NC_053024.1_cds_XP_048567703.1_14633__4572": "unresolved"
        }

    for leaf in t:
        if leaf.name in manual_candidates:
            leaf.props["two_odd_id"] = "candidate"
    _, ingroup_headers = split_seqs_by_2ODD_membership(seq_to_2ODD_id=seq_to_2ODD_id, tree=t)

    ingroup_headers = ingroup_headers - set(manual_candidates.keys())
    clusters = get_clusters(t)

    
    expected = {'0': {'two_odd_id': 'candidates_only',
                'percentage_of_ingroup_2ODD': None,
                'num_ingroup_2ODD': 0,
                'num_candidates': 1,
                'plant_groups': ['Dicots']},
                '1': {'two_odd_id': '2ODD34',
                'percentage_of_ingroup_2ODD': 0.026,
                'num_ingroup_2ODD': 1,
                'num_candidates': 0,
                'plant_groups': ['Basal Angiosperms']},
                '2': {'two_odd_id': '2ODD34',
                'percentage_of_ingroup_2ODD': 0.237,
                'num_ingroup_2ODD': 9,
                'num_candidates': 0,
                'plant_groups': ['Basal Angiosperms','Dicots']},
                '3': {'two_odd_id': '2ODD34',
                'percentage_of_ingroup_2ODD': 0.026,
                'num_ingroup_2ODD': 1,
                'num_candidates': 0,
                'plant_groups': ['Basal Angiosperms']},
                '4': {'two_odd_id': '2ODD34',
                'percentage_of_ingroup_2ODD': 0.026,
                'num_ingroup_2ODD': 1,
                'num_candidates': 0,
                'plant_groups': ['Dicots']},
                '5': {'two_odd_id': '2ODD34',
                'percentage_of_ingroup_2ODD': 0.105,
                'num_ingroup_2ODD': 4,
                'num_candidates': 0,
                'plant_groups': ['Basal Angiosperms','Dicots']},
                '6': {'two_odd_id': '2ODD34',
                'percentage_of_ingroup_2ODD': 0.579,
                'num_ingroup_2ODD': 22,
                'num_candidates': 9,
                'plant_groups': ['Dicots']},
                '7': {'two_odd_id': 'minor_2ODD_cluster',
                'percentage_of_ingroup_2ODD': None,
                'num_ingroup_2ODD': 3,
                'num_candidates': 1,
                'plant_groups': ['Dicots']},
                '8': {'two_odd_id': '2ODD38',
                'percentage_of_ingroup_2ODD': 0.029,
                'num_ingroup_2ODD': 3,
                'num_candidates': 0,
                'plant_groups': ['Dicots']},
                '9': {'two_odd_id': '2ODD38',
                'percentage_of_ingroup_2ODD': 0.039,
                'num_ingroup_2ODD': 4,
                'num_candidates': 0,
                'plant_groups': ['Dicots']},
                '10': {'two_odd_id': '2ODD38',
                'percentage_of_ingroup_2ODD': 0.01,
                'num_ingroup_2ODD': 1,
                'num_candidates': 0,
                'plant_groups': ['Dicots']},
                '11': {'two_odd_id': 'minor_2ODD_cluster',
                'percentage_of_ingroup_2ODD': None,
                'num_ingroup_2ODD': 1,
                'num_candidates': 0,
                'plant_groups': ['Dicots']},
                '12': {'two_odd_id': '2ODD37',
                'percentage_of_ingroup_2ODD': 0.067,
                'num_ingroup_2ODD': 2,
                'num_candidates': 0,
                'plant_groups': ['Dicots']},
                '13': {'two_odd_id': '2ODD37',
                'percentage_of_ingroup_2ODD': 0.033,
                'num_ingroup_2ODD': 1,
                'num_candidates': 0,
                'plant_groups': ['Dicots']},
                '14': {'two_odd_id': '2ODD37',
                'percentage_of_ingroup_2ODD': 0.567,
                'num_ingroup_2ODD': 17,
                'num_candidates': 0,
                'plant_groups': ['Dicots']},
                '15': {'two_odd_id': 'minor_2ODD_cluster',
                'percentage_of_ingroup_2ODD': None,
                'num_ingroup_2ODD': 9,
                'num_candidates': 0,
                'plant_groups': ['Dicots']},
                '16': {'two_odd_id': '2ODD37',
                'percentage_of_ingroup_2ODD': 0.033,
                'num_ingroup_2ODD': 1,
                'num_candidates': 0,
                'plant_groups': ['Dicots']},
                '17': {'two_odd_id': 'minor_2ODD_cluster',
                'percentage_of_ingroup_2ODD': None,
                'num_ingroup_2ODD': 1,
                'num_candidates': 0,
                'plant_groups': ['Dicots']},
                '18': {'two_odd_id': 'minor_2ODD_cluster',
                'percentage_of_ingroup_2ODD': None,
                'num_ingroup_2ODD': 2,
                'num_candidates': 0,
                'plant_groups': ['Dicots']},
                '19': {'two_odd_id': '2ODD37',
                'percentage_of_ingroup_2ODD': 0.3,
                'num_ingroup_2ODD': 9,
                'num_candidates': 0,
                'plant_groups': ['Dicots']},
                '20': {'two_odd_id': '2ODD35',
                'percentage_of_ingroup_2ODD': 0.326,
                'num_ingroup_2ODD': 15,
                'num_candidates': 6,
                'plant_groups': ['Dicots']},
                '21': {'two_odd_id': 'minor_2ODD_cluster',
                'percentage_of_ingroup_2ODD': None,
                'num_ingroup_2ODD': 4,
                'num_candidates': 0,
                'plant_groups': ['Dicots']},
                '22': {'two_odd_id': 'minor_2ODD_cluster',
                'percentage_of_ingroup_2ODD': None,
                'num_ingroup_2ODD': 14,
                'num_candidates': 0,
                'plant_groups': ['Dicots']},
                '23': {'two_odd_id': 'minor_2ODD_cluster',
                'percentage_of_ingroup_2ODD': None,
                'num_ingroup_2ODD': 6,
                'num_candidates': 0,
                'plant_groups': ['Dicots']},
                '24': {'two_odd_id': '2ODD36',
                'percentage_of_ingroup_2ODD': 0.628,
                'num_ingroup_2ODD': 49,
                'num_candidates': 1,
                'plant_groups': ['Dicots']},
                '25': {'two_odd_id': '2ODD36',
                'percentage_of_ingroup_2ODD': 0.333,
                'num_ingroup_2ODD': 26,
                'num_candidates': 0,
                'plant_groups': ['Dicots']},
                '26': {'two_odd_id': '2ODD36',
                'percentage_of_ingroup_2ODD': 0.038,
                'num_ingroup_2ODD': 3,
                'num_candidates': 0,
                'plant_groups': ['Dicots']},
                '27': {'two_odd_id': '2ODD35',
                'percentage_of_ingroup_2ODD': 0.674,
                'num_ingroup_2ODD': 31,
                'num_candidates': 0,
                'plant_groups': ['Dicots']},
                '28': {'two_odd_id': '2ODD38',
                'percentage_of_ingroup_2ODD': 0.922,
                'num_ingroup_2ODD': 95,
                'num_candidates': 1,
                'plant_groups': ['Dicots', 'Monocots']}}

    result = cluster_meta_info(clusters, major_minor_2ODDs_dict,ingroup_headers, manual_candidates)
    assert result == expected

def test_get_candidate_to_char_baits_df():
    manual_candidates = {
        "lcl_NC_084852.1_cds_XP_061960601.1_3168__3691" : "2ODD34",
        'VIT_204s0008g04920.2__29760' : '2ODD34',
        'Ptrif.0001s0318.1__37690' : '2ODD34',
        'Atru_chr7_2342__47965' : '2ODD34',
        'lcl_NW_019168159.1_cds_XP_022728803.1_47147__66656' : '2ODD34',
        'lcl_NC_045134.1_cds_XP_031407529.1_33286__22663' : '2ODD34',
        'GWHTACBH000860__413952' : '2ODD34',
        'Potri.006G248000.1__3694' : '2ODD34',
        'lcl_NC_065570.1_cds_XP_050209042.1_2787__3986' : '2ODD34',
        'Lus10008097_PACid-23169790__4006' : '2ODD34', 
        'Acc09929.1__3625': "minor_2ODD_cluster",
        'lcl_NC_031989.1_cds_XP_019252277.1_2282__49451': "2ODD35",
        'Kaladp0011s0492.1.v1.1__63787': "2ODD35",
        'lcl_OX459121.1_cds_CAI9101307.1_14231__43536': "2ODD35",
        'rna-gnl_WGS-JAMLDZ_EVM0034906.1__4058': "2ODD35",
        'lcl_NC_039901.1_cds_XP_027112943.1_15453__13443': "2ODD35",
        'rna-gnl_WGS-JAMLDZ_EVM0007978.1__4058': "2ODD35", 
        'Bol038153__3712': "2ODD36",
        "lcl_NC_053024.1_cds_XP_048567703.1_14633__4572": "unresolved"
        }

    for leaf in t:
        if leaf.name in manual_candidates:
            leaf.props["two_odd_id"] = "candidate"
    _, ingroup_headers = split_seqs_by_2ODD_membership(seq_to_2ODD_id=seq_to_2ODD_id, tree=t) 
    seq_id_to_idx_dict = seq_to_cluster_idx(get_clusters(t))
    expected_df = pd.DataFrame({'candidate': {0: 'lcl_NC_084852.1_cds_XP_061960601.1_3168__3691',
                1: 'VIT_204s0008g04920.2__29760',
                2: 'Ptrif.0001s0318.1__37690',
                3: 'Atru_chr7_2342__47965',
                4: 'lcl_NW_019168159.1_cds_XP_022728803.1_47147__66656',
                5: 'lcl_NC_045134.1_cds_XP_031407529.1_33286__22663',
                6: 'GWHTACBH000860__413952',
                7: 'Potri.006G248000.1__3694',
                8: 'lcl_NC_065570.1_cds_XP_050209042.1_2787__3986',
                9: 'Lus10008097_PACid-23169790__4006',
                10: 'Acc09929.1__3625',
                11: 'lcl_NC_031989.1_cds_XP_019252277.1_2282__49451',
                12: 'Kaladp0011s0492.1.v1.1__63787',
                13: 'lcl_OX459121.1_cds_CAI9101307.1_14231__43536',
                14: 'rna-gnl_WGS-JAMLDZ_EVM0034906.1__4058',
                15: 'lcl_NC_039901.1_cds_XP_027112943.1_15453__13443',
                16: 'rna-gnl_WGS-JAMLDZ_EVM0007978.1__4058',
                17: 'Bol038153__3712',
                18: 'lcl_NC_053024.1_cds_XP_048567703.1_14633__4572'},
                'cluster_idx': {0: 0,
                1: 6,
                2: 6,
                3: 6,
                4: 6,
                5: 6,
                6: 6,
                7: 6,
                8: 6,
                9: 6,
                10: 7,
                11: 20,
                12: 20,
                13: 20,
                14: 20,
                15: 20,
                16: 20,
                17: 24,
                18: 28},
                'closest_char_bait': {0: 'AF417859__AOP3__glucosinolate_biosynthesis__3702',
                1: 'AF417859__AOP3__glucosinolate_biosynthesis__3702',
                2: 'AF417859__AOP3__glucosinolate_biosynthesis__3702',
                3: 'AF417859__AOP3__glucosinolate_biosynthesis__3702',
                4: 'AF417859__AOP3__glucosinolate_biosynthesis__3702',
                5: 'AF417859__AOP3__glucosinolate_biosynthesis__3702',
                6: 'AF417859__AOP3__glucosinolate_biosynthesis__3702',
                7: 'AF417859__AOP3__glucosinolate_biosynthesis__3702',
                8: 'AF417859__AOP3__glucosinolate_biosynthesis__3702',
                9: 'AF417859__AOP3__glucosinolate_biosynthesis__3702',
                10: 'AF417859__AOP3__glucosinolate_biosynthesis__3702',
                11: 'AF417859__AOP3__glucosinolate_biosynthesis__3702',
                12: 'AF417859__AOP3__glucosinolate_biosynthesis__3702',
                13: 'AF417859__AOP3__glucosinolate_biosynthesis__3702',
                14: 'AF417859__AOP3__glucosinolate_biosynthesis__3702',
                15: 'AF417859__AOP3__glucosinolate_biosynthesis__3702',
                16: 'AF417859__AOP3__glucosinolate_biosynthesis__3702',
                17: 'AF417859__AOP3__glucosinolate_biosynthesis__3702',
                18: 'AF417859__AOP3__glucosinolate_biosynthesis__3702'},
                'closest_char_bait_dist': {0: 2.3441771,
                1: 2.2212721999999996,
                2: 2.4032200999999995,
                3: 2.4297531,
                4: 2.2819316999999995,
                5: 2.4276866999999998,
                6: 2.296112,
                7: 2.2915937,
                8: 2.4520146999999994,
                9: 2.6844072,
                10: 2.7427177,
                11: 2.1061180999999998,
                12: 2.1150531,
                13: 2.2243331,
                14: 2.2426211,
                15: 2.2775660999999996,
                16: 2.2996330999999994,
                17: 2.5819186,
                18: 3.9759043999999997},
                'second_closest_char_bait': {0: 'AF417858__AOP2__glucosinolate_biosynthesis__3702',
                1: 'AF417858__AOP2__glucosinolate_biosynthesis__3702',
                2: 'AF417858__AOP2__glucosinolate_biosynthesis__3702',
                3: 'AF417858__AOP2__glucosinolate_biosynthesis__3702',
                4: 'AF417858__AOP2__glucosinolate_biosynthesis__3702',
                5: 'AF417858__AOP2__glucosinolate_biosynthesis__3702',
                6: 'AF417858__AOP2__glucosinolate_biosynthesis__3702',
                7: 'AF417858__AOP2__glucosinolate_biosynthesis__3702',
                8: 'AF417858__AOP2__glucosinolate_biosynthesis__3702',
                9: 'AF417858__AOP2__glucosinolate_biosynthesis__3702',
                10: 'AF417858__AOP2__glucosinolate_biosynthesis__3702',
                11: 'AF417858__AOP2__glucosinolate_biosynthesis__3702',
                12: 'AF417858__AOP2__glucosinolate_biosynthesis__3702',
                13: 'AF417858__AOP2__glucosinolate_biosynthesis__3702',
                14: 'AF417858__AOP2__glucosinolate_biosynthesis__3702',
                15: 'AF417858__AOP2__glucosinolate_biosynthesis__3702',
                16: 'AF417858__AOP2__glucosinolate_biosynthesis__3702',
                17: 'AF417858__AOP2__glucosinolate_biosynthesis__3702',
                18: 'AF417858__AOP2__glucosinolate_biosynthesis__3702'},
                'second_closest_char_bait_dist': {0: 2.3673281,
                1: 2.2444232,
                2: 2.4263711,
                3: 2.4529041,
                4: 2.3050827,
                5: 2.4508377,
                6: 2.319263,
                7: 2.3147447,
                8: 2.4751657,
                9: 2.7075582,
                10: 2.7658687,
                11: 2.1292690999999997,
                12: 2.1382041,
                13: 2.2474841,
                14: 2.2657721,
                15: 2.3007171,
                16: 2.3227841,
                17: 2.6050695999999998,
                18: 3.9990553999999996
                }
            }
    )

    result = get_candidate_to_char_baits_df(manual_candidates.keys(), ingroup_headers, dist_dict, seq_id_to_idx_dict)
    assert result.equals(expected_df)

