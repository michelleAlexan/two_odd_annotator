#%%
import json

import pytest
from pathlib import Path
import random
from ete4 import PhyloTree, Tree

from bio_tools.viz.tree import (
    explore_tree_plant_groups,
    load_treecluster_assignments, 
    explore_tree_cluster_clades
    )


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
    get_landscape, 
    resolve_candidates_in_landscape,
    two_odd_id_to_landscape_indices, 
    seq_id_to_landscape_idx,
    landscape_meta_info)

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


#%%
RESULTS_DIR = Path(__file__).parents[1] /  "results" 

config = load_config(Path(__file__).parents[2] / DEFAULT_CONFIG_PATH)
config["annotate"]["ingroup"] = Path(__file__).parents[2] / "data" / "2ODDs" / "characterized_2ODDs.fasta"

major_minor_2ODDs_path = Path(__file__).parents[2] / "data" / "2ODDs" / "major_minor_2ODD_ids_manual.json"
major_minor_2ODDs_dict = json.load(open(major_minor_2ODDs_path))

seq_to_2ODD_id = reverse_major_minor_2ODD_dict(major_minor_2ODD_dict=major_minor_2ODDs_dict)




def test_create_annotation_fasta():
    if (RESULTS_DIR / ANNOTATION_FASTA).exists():
        (RESULTS_DIR / ANNOTATION_FASTA).unlink()   

    
    create_annotation_fasta(
        results_dir=RESULTS_DIR,
        ingroup_2ODD_fasta=config["annotate"]["ingroup"],
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

    assert t_char_2ODDs["PV023584__F3H__flavonoid_pathway__981085"].props["two_odd_id"] == "2ODD14"
    assert t_char_2ODDs["sp|Q96323.1_ANS_Arabidopsis_thaliana__3702"].props["two_odd_id"] == "candidate"

def test_assign_plant_group_props():
    assign_plant_group_props(tree=t_char_2ODDs)

    assert t_char_2ODDs["At4g10500__S3H__salicylic_acid_metabolism__3702"].props["plant_group"] == "Dicots"
    assert t_char_2ODDs["Os04g49210__S5H__salicylic_acid_metabolism__4530"].props["plant_group"] == "Monocots"


#%% test landscape clustering functions

test_anno_tree_path = Path(__file__).parents[1] / "data" / "test_anno_tree.nwk"
t = PhyloTree(open(test_anno_tree_path), sp_naming_function=lambda name: name.split('__')[-1])
tax2names, tax2lineages, tax2rank = t.annotate_ncbi_taxa(taxid_attr='species')

candidate_headers, ingroup_headers = split_seqs_by_2ODD_membership(seq_to_2ODD_id=seq_to_2ODD_id, tree=t)


assign_2ODD_props(
    tree=t, seq_to_2ODD_id=seq_to_2ODD_id, 
    candidate_headers=candidate_headers
)
assign_plant_group_props(tree=t)
assign_plant_group_props(tree=t)

# explore_tree_cluster_clades(t)
#%% TEST LANDSCAPE CLUSTERING

expected_landscape_without_candidates = [
    [
        t["lcl_NC_084852.1_cds_XP_061960601.1_3168__3691"] # 2ODD41
    ], 
    [
        t["rna-XM_031626300.2__210225"], # 2ODD37
        t["CKAN_00595100__337451"], 
        t["Aristolochia_fimbriatalcl_CM034078.1_cds_KAG9459387.1_3895__158543"], 
        t["Gobar.D13G239800.1__3634"],
        t["Lj8A611G43.1__105884"], 
        t["FvH4_2g29990.t1__57918"],
        t['Casgl23S04257__3522'],
        t['CiLak.07G220400.1__32201'],
        t['lcl_NC_065570.1_cds_XP_050209022.1_2785__3986'],
        t['Lus10013130_PACid-23164132__4006'],
        t['evm_27.model.AmTr_v1.0_scaffold00174.4__13333'],
        t['RZC71586__3469'],
        t['lcl_NW_010729080.1_cds_XP_010249917.1_5454__4432'],
        t['rna-XM_042628350.1__60698'],
        t['lcl_CM056811.1_cds_KAJ8637692.1_11945__3435'],
        t['Aristolochia_fimbriatalcl_CM034078.1_cds_KAG9459386.1_3894__158543'],
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
        t["Acc09929.1__3625"], # 2ODD37
    ], 
    [
        t["pveT_jg33784.t1__170927"] # minor 2ODD cluster
    ], 
    [
        t['Eucgr.C03724.1__71139'], # 2ODD41
        t['lcl_NC_045130.1_cds_XP_031390300.1_19564__22663'],
        t['Eucgr.C04147.1__71139'],
        t['Acc23206.1__3625'],
        t['lcl_CM027379.1_cds_KAF9612985.1_11655__261450'],
        t['RZC71605__3469'],
        t['rna-XM_042647017.1__60698'],
        t['lcl_CM027379.1_cds_KAF9612671.1_9820__261450']
    ],
    [
        t['Lj7A541T77.1__105884'] # minor 2ODD cluster
    ],
    [
        t['Pg_S7309.1__4054'],  # 2ODD40
        t['DCAR_015859__4039'],
        t['CsatW809539.1__3659'], 
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
        t['pveT_jg9723.t1__170927'] # 2ODD40
    ],
    [
        t["Cc06_g15940.1__49390"], # minor 2ODD cluster
        t["Gobar.D05G137100.1__3634"],
        t["lcl_OX459118.1_cds_CAI9089985.1_2909__43536"],
    ],
    [
        t['TRINITY_DN197356_c3_g1_i1__4102'], # 2ODD40
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
    ],
    [
        t["Kaladp0011s0492.1.v1.1__63787"],
        t['lcl_OX459121.1_cds_CAI9101307.1_14231__43536'], # 2ODD40
        t['rna-gnl_WGS-JAMLDZ_EVM0034906.1__4058'], 
        t['lcl_NC_039901.1_cds_XP_027112943.1_15453__13443'],
        t["rna-gnl_WGS-JAMLDZ_EVM0007978.1__4058"]
    ], 
    [
        t['lcl_NC_080155.1_cds_XP_057473373.1_43266__165200'], # 2ODD38
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
        t['lcl_NC_084796.1_cds_XP_062119410.1_20878__3486'],
        t['Umino11360.1__262084'],
        t['FvH4_3g36530.t1__57918'],
        t['lcl_NC_083380.1_cds_XP_015881543.2_5733__326968'],
        t['Kaladp0007s0029.1.v1.1__63787'],
        t['Blora12794.1__200023']
    ], 
    [
        t['VIT_209s0002g05340.1__29760'], # 2ODD39
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
        t['CiLak.02G150500.1__32201'],
        t['FSB014894201__28930'],
        t['lcl_NC_044908.1_cds_XP_030968569.1_22136__97700']
 ], 
 [
        t['lcl_PKPP01002761.1_cds_PWA73280.1_26522__35608'], # 2ODD38
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
        t['lcl_CM008454.1_cds_PHT32940.1_29277__33114'], # 2ODD41
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
        t['strangu_020525-RA__38733']
    ], 
    [
        t['lcl_NC_053024.1_cds_XP_048567703.1_14633__4572']  # minor 2ODD cluster
    ]
]

def test_get_landscape_without_candidates():
    result = get_landscape(t)
    assert result == expected_landscape_without_candidates



#%%
# take the test tree with all input sequences, manually select candidates, and create a dictionary with manually assigned candidates as keys that map to the expected output of the annotation service.
# map candidate name to expected cluster id resulted from the annotation service.
# The expected cluster id IS NOT the cluster id of the actual input sequence 

expected_landscape_with_candidates = [
    [
        t["lcl_NC_084852.1_cds_XP_061960601.1_3168__3691"] # candidate
    ], 
    [
        t["rna-XM_031626300.2__210225"], # 2ODD37
        t["CKAN_00595100__337451"], 
        t["Aristolochia_fimbriatalcl_CM034078.1_cds_KAG9459387.1_3895__158543"], 
        t["Gobar.D13G239800.1__3634"],
        t["Lj8A611G43.1__105884"], 
        t["FvH4_2g29990.t1__57918"],
        t['Casgl23S04257__3522'],
        t['CiLak.07G220400.1__32201'],
        t['lcl_NC_065570.1_cds_XP_050209022.1_2785__3986'],
        t['Lus10013130_PACid-23164132__4006'],
        t['evm_27.model.AmTr_v1.0_scaffold00174.4__13333'],
        t['RZC71586__3469'],
        t['lcl_NW_010729080.1_cds_XP_010249917.1_5454__4432'],
        t['rna-XM_042628350.1__60698'],
        t['lcl_CM056811.1_cds_KAJ8637692.1_11945__3435'],
        t['Aristolochia_fimbriatalcl_CM034078.1_cds_KAG9459386.1_3894__158543'],
        t['Kaladp1006s0012.1.v1.1__63787']

    ],
    [        
        t['VIT_204s0008g04920.2__29760'], # candidate 
        t['Ptrif.0001s0318.1__37690'],
        t['Atru_chr7_2342__47965'],
        t['lcl_NW_019168159.1_cds_XP_022728803.1_47147__66656'],
        t['lcl_NC_045134.1_cds_XP_031407529.1_33286__22663'],
        t['GWHTACBH000860__413952'],
        t['Potri.006G248000.1__3694'],
        t['lcl_NC_065570.1_cds_XP_050209042.1_2787__3986'],
        t['Lus10008097_PACid-23169790__4006']    
    ],
    [

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
        t["Acc09929.1__3625"], # candidate
    ], 
    [
        t["pveT_jg33784.t1__170927"] # minor 2ODD cluster
    ], 
    [
        t['Eucgr.C03724.1__71139'], # 2ODD41
        t['lcl_NC_045130.1_cds_XP_031390300.1_19564__22663'],
        t['Eucgr.C04147.1__71139'],
        t['Acc23206.1__3625'],
        t['lcl_CM027379.1_cds_KAF9612985.1_11655__261450'],
        t['RZC71605__3469'],
        t['rna-XM_042647017.1__60698'],
        t['lcl_CM027379.1_cds_KAF9612671.1_9820__261450']
    ],
    [
        t['Lj7A541T77.1__105884'] # minor 2ODD cluster
    ],
    [
        t['Pg_S7309.1__4054'],  # 2ODD40
        t['DCAR_015859__4039'],
        t['CsatW809539.1__3659'], 
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
        t['pveT_jg9723.t1__170927'] # 2ODD40
    ],
    [
        t["Cc06_g15940.1__49390"], # minor 2ODD cluster
        t["Gobar.D05G137100.1__3634"],
        t["lcl_OX459118.1_cds_CAI9089985.1_2909__43536"],
    ],
    [
        t['TRINITY_DN197356_c3_g1_i1__4102'], # 2ODD40
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
        t["lcl_NC_031989.1_cds_XP_019252277.1_2282__49451"], # candidate
        t["Kaladp0011s0492.1.v1.1__63787"],
        t['lcl_OX459121.1_cds_CAI9101307.1_14231__43536'], 
        t['rna-gnl_WGS-JAMLDZ_EVM0034906.1__4058'], 
        t['lcl_NC_039901.1_cds_XP_027112943.1_15453__13443'],
        t["rna-gnl_WGS-JAMLDZ_EVM0007978.1__4058"]
    ], 
    [
        t['lcl_NC_080155.1_cds_XP_057473373.1_43266__165200'], # 2ODD38
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
        t['lcl_NC_084796.1_cds_XP_062119410.1_20878__3486'],
        t['Umino11360.1__262084'],
        t['FvH4_3g36530.t1__57918'],
        t['lcl_NC_083380.1_cds_XP_015881543.2_5733__326968'],
        t['Kaladp0007s0029.1.v1.1__63787'],
        t['Blora12794.1__200023']
    ], 
    [
        t['VIT_209s0002g05340.1__29760'], # 2ODD39
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
    ],
    [
        t['Bol038153__3712'], # candidate
    ],
    [
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
        t['CiLak.02G150500.1__32201'],
        t['FSB014894201__28930'],
        t['lcl_NC_044908.1_cds_XP_030968569.1_22136__97700']
    ], 
    [
        t['lcl_PKPP01002761.1_cds_PWA73280.1_26522__35608'], # 2ODD38
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
        t['lcl_CM008454.1_cds_PHT32940.1_29277__33114'], # 2ODD41
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
        t['strangu_020525-RA__38733']
    ], 
    [
        t['lcl_NC_053024.1_cds_XP_048567703.1_14633__4572']  # candidate
    ]
]

def test_landscape_with_candidates():
    manual_candidates = {
    "lcl_NC_084852.1_cds_XP_061960601.1_3168__3691" : "2ODD37",
    'VIT_204s0008g04920.2__29760' : '2ODD37',
    'Ptrif.0001s0318.1__37690' : '2ODD37',
    'Atru_chr7_2342__47965' : '2ODD37',
    'lcl_NW_019168159.1_cds_XP_022728803.1_47147__66656' : '2ODD37',
    'lcl_NC_045134.1_cds_XP_031407529.1_33286__22663' : '2ODD37',
    'GWHTACBH000860__413952' : '2ODD37',
    'Potri.006G248000.1__3694' : '2ODD37',
    'lcl_NC_065570.1_cds_XP_050209042.1_2787__3986' : '2ODD37',
    'Lus10008097_PACid-23169790__4006' : '2ODD37', 
    'Acc09929.1__3625': "minor_2ODD_cluster",
    'lcl_NC_031989.1_cds_XP_019252277.1_2282__49451': "2ODD38",
    'Kaladp0011s0492.1.v1.1__63787': "2ODD38",
    'lcl_OX459121.1_cds_CAI9101307.1_14231__43536': "2ODD38",
    'rna-gnl_WGS-JAMLDZ_EVM0034906.1__4058': "2ODD38",
    'lcl_NC_039901.1_cds_XP_027112943.1_15453__13443': "2ODD38",
    'rna-gnl_WGS-JAMLDZ_EVM0007978.1__4058': "2ODD38", 
    'Bol038153__3712': "2ODD39",
    "lcl_NC_053024.1_cds_XP_048567703.1_14633__4572": "unresolved"
    }

    for leaf in t:
        if leaf.name in manual_candidates:
            leaf.props["two_odd_id"] = "candidate"

    result = get_landscape(t)
    assert result == expected_landscape_with_candidates




#%%  test the resolved landscape 
dist_dict = build_distance_lookup(t)

#get median distance between all pairs of leaves in the tree
all_distances = []
for leaf1 in t.leaves():
    for leaf2 in t.leaves():
        if leaf1 != leaf2:
            all_distances.append(dist_dict[leaf1.name][leaf2.name])
median_distance = sorted(all_distances)[len(all_distances) // 2]


expected_resolved_landscape_with_candidates = [
    [
        t["lcl_NC_084852.1_cds_XP_061960601.1_3168__3691"],
        t["rna-XM_031626300.2__210225"], # 2ODD37
        t["CKAN_00595100__337451"], 
        t["Aristolochia_fimbriatalcl_CM034078.1_cds_KAG9459387.1_3895__158543"], 
        t["Gobar.D13G239800.1__3634"],
        t["Lj8A611G43.1__105884"], 
        t["FvH4_2g29990.t1__57918"],
        t['Casgl23S04257__3522'],
        t['CiLak.07G220400.1__32201'],
        t['lcl_NC_065570.1_cds_XP_050209022.1_2785__3986'],
        t['Lus10013130_PACid-23164132__4006'],
        t['evm_27.model.AmTr_v1.0_scaffold00174.4__13333'],
        t['RZC71586__3469'],
        t['lcl_NW_010729080.1_cds_XP_010249917.1_5454__4432'],
        t['rna-XM_042628350.1__60698'],
        t['lcl_CM056811.1_cds_KAJ8637692.1_11945__3435'],
        t['Aristolochia_fimbriatalcl_CM034078.1_cds_KAG9459386.1_3894__158543'],
        t['Kaladp1006s0012.1.v1.1__63787'],
        t['VIT_204s0008g04920.2__29760'], # candidate 
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
        t["Acc09929.1__3625"],
        t["pveT_jg33784.t1__170927"] # minor 2ODD cluster
    ], 
    [
        t['Eucgr.C03724.1__71139'], # 2ODD41
        t['lcl_NC_045130.1_cds_XP_031390300.1_19564__22663'],
        t['Eucgr.C04147.1__71139'],
        t['Acc23206.1__3625'],
        t['lcl_CM027379.1_cds_KAF9612985.1_11655__261450'],
        t['RZC71605__3469'],
        t['rna-XM_042647017.1__60698'],
        t['lcl_CM027379.1_cds_KAF9612671.1_9820__261450']
    ],
    [
        t['Lj7A541T77.1__105884'] # minor 2ODD cluster
    ],
    [
        t['Pg_S7309.1__4054'],  # 2ODD40
        t['DCAR_015859__4039'],
        t['CsatW809539.1__3659'], 
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
        t['pveT_jg9723.t1__170927'] # 2ODD40
    ],
    [
        t["Cc06_g15940.1__49390"], # minor 2ODD cluster
        t["Gobar.D05G137100.1__3634"],
        t["lcl_OX459118.1_cds_CAI9089985.1_2909__43536"],
    ],
    [
        t['TRINITY_DN197356_c3_g1_i1__4102'], # 2ODD40
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
        t["lcl_NC_031989.1_cds_XP_019252277.1_2282__49451"], # 2ODD38
        t["Kaladp0011s0492.1.v1.1__63787"],
        t['lcl_OX459121.1_cds_CAI9101307.1_14231__43536'], 
        t['rna-gnl_WGS-JAMLDZ_EVM0034906.1__4058'], 
        t['lcl_NC_039901.1_cds_XP_027112943.1_15453__13443'],
        t["rna-gnl_WGS-JAMLDZ_EVM0007978.1__4058"],
        t['lcl_NC_080155.1_cds_XP_057473373.1_43266__165200'], 
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
        t['lcl_NC_084796.1_cds_XP_062119410.1_20878__3486'],
        t['Umino11360.1__262084'],
        t['FvH4_3g36530.t1__57918'],
        t['lcl_NC_083380.1_cds_XP_015881543.2_5733__326968'],
        t['Kaladp0007s0029.1.v1.1__63787'],
        t['Blora12794.1__200023']
    ], 
    [
        t['VIT_209s0002g05340.1__29760'], # 2ODD39
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
        t['CiLak.02G150500.1__32201'],
        t['FSB014894201__28930'],
        t['lcl_NC_044908.1_cds_XP_030968569.1_22136__97700']
 ], 
 [
        t['lcl_PKPP01002761.1_cds_PWA73280.1_26522__35608'], # 2ODD38
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
        t['lcl_CM008454.1_cds_PHT32940.1_29277__33114'], # 2ODD41
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
        t['strangu_020525-RA__38733']
    ], 
    [
        t['lcl_NC_053024.1_cds_XP_048567703.1_14633__4572']  # unresolved
    ]
]

def test_resolved_landscape_with_candidates():
    landscape_unresolved = get_landscape(t)
    resolved_landscape = resolve_candidates_in_landscape(landscape=landscape_unresolved, 
                                                         dist_dict=dist_dict, 
                                                         threshold=2)

    assert resolved_landscape == expected_resolved_landscape_with_candidates


#%% 


def test_two_odd_id_to_landscape_indices():
    manual_candidates = {
        "lcl_NC_084852.1_cds_XP_061960601.1_3168__3691" : "2ODD37",
        'VIT_204s0008g04920.2__29760' : '2ODD37',
        'Ptrif.0001s0318.1__37690' : '2ODD37',
        'Atru_chr7_2342__47965' : '2ODD37',
        'lcl_NW_019168159.1_cds_XP_022728803.1_47147__66656' : '2ODD37',
        'lcl_NC_045134.1_cds_XP_031407529.1_33286__22663' : '2ODD37',
        'GWHTACBH000860__413952' : '2ODD37',
        'Potri.006G248000.1__3694' : '2ODD37',
        'lcl_NC_065570.1_cds_XP_050209042.1_2787__3986' : '2ODD37',
        'Lus10008097_PACid-23169790__4006' : '2ODD37', 
        'Acc09929.1__3625': "minor_2ODD_cluster",
        'lcl_NC_031989.1_cds_XP_019252277.1_2282__49451': "2ODD38",
        'Kaladp0011s0492.1.v1.1__63787': "2ODD38",
        'lcl_OX459121.1_cds_CAI9101307.1_14231__43536': "2ODD38",
        'rna-gnl_WGS-JAMLDZ_EVM0034906.1__4058': "2ODD38",
        'lcl_NC_039901.1_cds_XP_027112943.1_15453__13443': "2ODD38",
        'rna-gnl_WGS-JAMLDZ_EVM0007978.1__4058': "2ODD38", 
        'Bol038153__3712': "2ODD39",
        "lcl_NC_053024.1_cds_XP_048567703.1_14633__4572": "unresolved"
        }

    for leaf in t:
        if leaf.name in manual_candidates:
            leaf.props["two_odd_id"] = "candidate"

    landscape_unresolved = get_landscape(t)
    resolved_landscape = resolve_candidates_in_landscape(landscape=landscape_unresolved, 
                                                            dist_dict=dist_dict, 
                                                            threshold=2)

    expected = {'2ODD37': [0],
                'minor_2ODD_cluster': [1, 3, 5, 7, 10],
                '2ODD41': [2, 13],
                '2ODD40': [4, 6, 8],
                '2ODD38': [9, 12],
                '2ODD39': [11],
                'unresolved': [14]}

    result = two_odd_id_to_landscape_indices(resolved_landscape)
    assert result == expected


def test_seq_id_to_landscape_idx():

    manual_candidates = {
        "lcl_NC_084852.1_cds_XP_061960601.1_3168__3691" : "2ODD37",
        'VIT_204s0008g04920.2__29760' : '2ODD37',
        'Ptrif.0001s0318.1__37690' : '2ODD37',
        'Atru_chr7_2342__47965' : '2ODD37',
        'lcl_NW_019168159.1_cds_XP_022728803.1_47147__66656' : '2ODD37',
        'lcl_NC_045134.1_cds_XP_031407529.1_33286__22663' : '2ODD37',
        'GWHTACBH000860__413952' : '2ODD37',
        'Potri.006G248000.1__3694' : '2ODD37',
        'lcl_NC_065570.1_cds_XP_050209042.1_2787__3986' : '2ODD37',
        'Lus10008097_PACid-23169790__4006' : '2ODD37', 
        'Acc09929.1__3625': "minor_2ODD_cluster",
        'lcl_NC_031989.1_cds_XP_019252277.1_2282__49451': "2ODD38",
        'Kaladp0011s0492.1.v1.1__63787': "2ODD38",
        'lcl_OX459121.1_cds_CAI9101307.1_14231__43536': "2ODD38",
        'rna-gnl_WGS-JAMLDZ_EVM0034906.1__4058': "2ODD38",
        'lcl_NC_039901.1_cds_XP_027112943.1_15453__13443': "2ODD38",
        'rna-gnl_WGS-JAMLDZ_EVM0007978.1__4058': "2ODD38", 
        'Bol038153__3712': "2ODD39",
        "lcl_NC_053024.1_cds_XP_048567703.1_14633__4572": "unresolved"
        }

    for leaf in t:
        if leaf.name in manual_candidates:
            leaf.props["two_odd_id"] = "candidate"

    landscape_unresolved = get_landscape(t)
    resolved_landscape = resolve_candidates_in_landscape(landscape=landscape_unresolved, 
                                                            dist_dict=dist_dict, 
                                                            threshold=2)
    

    expected = {'lcl_NC_084852.1_cds_XP_061960601.1_3168__3691': 0,
        'rna-XM_031626300.2__210225': 0,
        'CKAN_00595100__337451': 0,
        'Aristolochia_fimbriatalcl_CM034078.1_cds_KAG9459387.1_3895__158543': 0,
        'Gobar.D13G239800.1__3634': 0,
        'Lj8A611G43.1__105884': 0,
        'FvH4_2g29990.t1__57918': 0,
        'Casgl23S04257__3522': 0,
        'CiLak.07G220400.1__32201': 0,
        'lcl_NC_065570.1_cds_XP_050209022.1_2785__3986': 0,
        'Lus10013130_PACid-23164132__4006': 0,
        'evm_27.model.AmTr_v1.0_scaffold00174.4__13333': 0,
        'RZC71586__3469': 0,
        'lcl_NW_010729080.1_cds_XP_010249917.1_5454__4432': 0,
        'rna-XM_042628350.1__60698': 0,
        'lcl_CM056811.1_cds_KAJ8637692.1_11945__3435': 0,
        'Aristolochia_fimbriatalcl_CM034078.1_cds_KAG9459386.1_3894__158543': 0,
        'Kaladp1006s0012.1.v1.1__63787': 0,
        'VIT_204s0008g04920.2__29760': 0,
        'Ptrif.0001s0318.1__37690': 0,
        'Atru_chr7_2342__47965': 0,
        'lcl_NW_019168159.1_cds_XP_022728803.1_47147__66656': 0,
        'lcl_NC_045134.1_cds_XP_031407529.1_33286__22663': 0,
        'GWHTACBH000860__413952': 0,
        'Potri.006G248000.1__3694': 0,
        'lcl_NC_065570.1_cds_XP_050209042.1_2787__3986': 0,
        'Lus10008097_PACid-23169790__4006': 0,
        'Acc06446.1__3625': 0,
        'Acc16585.1__3625': 0,
        'pveT_jg24143.t1__170927': 0,
        'pveT_jg28701.t1__170927': 0,
        'lcl_NC_083383.1_cds_XP_048332159.2_17459__326968': 0,
        'Cnepa31456.1__79760': 0,
        'Umino40998.1__262084': 0,
        'FvH4_2g29991.t1__57918': 0,
        'Lj1g0022518.1__34305': 0,
        'Lj4g0022664.1__34305': 0,
        'Medtr1g011600.1__3880': 0,
        'Ca_10933__3827': 0,
        'Ler.1DRT.1g009880.1__41257': 0,
        'Ler.1DRT.3g073930.1__41257': 0,
        'Medtr3g108520.1__3880': 0,
        'Ca_23071__3827': 0,
        'Ca_01791__3827': 0,
        'Medtr4g021360.1__3880': 0,
        'Ca_05192__3827': 0,
        'Medtr4g021380.1__3880': 0,
        'lcl_OX451737.1_cds_CAI8600332.1_14923__3906': 0,
        'pveT_jg20218.t1__170927': 1,
        'pveT_jg20220.t1__170927': 1,
        'Acc09929.1__3625': 1,
        'pveT_jg33784.t1__170927': 1,
        'Eucgr.C03724.1__71139': 2,
        'lcl_NC_045130.1_cds_XP_031390300.1_19564__22663': 2,
        'Eucgr.C04147.1__71139': 2,
        'Acc23206.1__3625': 2,
        'lcl_CM027379.1_cds_KAF9612985.1_11655__261450': 2,
        'RZC71605__3469': 2,
        'rna-XM_042647017.1__60698': 2,
        'lcl_CM027379.1_cds_KAF9612671.1_9820__261450': 2,
        'Lj7A541T77.1__105884': 3,
        'Pg_S7309.1__4054': 4,
        'DCAR_015859__4039': 4,
        'CsatW809539.1__3659': 4,
        'lcl_NW_026775580.1_cds_XP_059633894.1_39517__4283': 4,
        'Dglom27812.1__34297': 4,
        'lcl_CM059866.1_cds_KAK0595785.1_1317__4024': 4,
        'Ptrif.0007s1329.1__37690': 4,
        'Potri.001G176500.1__3694': 4,
        'lcl_NC_065574.1_cds_XP_050234707.1_20560__3986': 4,
        'Umino20162.1__262084': 4,
        'lcl_NC_083379.1_cds_XP_015868728.2_4750__326968': 4,
        'lcl_NC_084800.1_cds_XP_062086870.1_35676__3486': 4,
        'lcl_NC_045131.1_cds_XP_031396687.1_21997__22663': 4,
        'Eucgr.J03022.1__71139': 4,
        'VIT_209s0002g05290.1__29760': 4,
        'lcl_NW_026775580.1_cds_XP_059636056.1_39514__4283': 4,
        'Acc08999.1__3625': 4,
        'Acc07679.1__3625': 4,
        'FSB010835101__28930': 4,
        'DCAR_015858__4039': 4,
        'Casgl236S24615__3522': 5,
        'CiLak.15G093300.1__32201': 5,
        'Qurub.10G093800.1__3512': 5,
        'lcl_CM025849.1_cds_KAB1226463.1_1752__262757': 5,
        'Tcord05348.2__703396': 5,
        'HU01G00826.1__176265': 5,
        'Sc13g0006470.01__3999': 5,
        'TRINITY_DN108901_c0_g1_i1__223224': 5,
        'CM008280.1.CM008280.1.g987.t1__62330': 5,
        'pveT_jg9723.t1__170927': 6,
        'Cc06_g15940.1__49390': 7,
        'Gobar.D05G137100.1__3634': 7,
        'lcl_OX459118.1_cds_CAI9089985.1_2909__43536': 7,
        'TRINITY_DN197356_c3_g1_i1__4102': 8,
        'lcl_CM061488.1_cds_KAK1430368.1_13046__13708': 8,
        'Oeu000043.1__4146': 8,
        'lcl_NW_026137546.1_cds_XP_051123007.1_5823__175694': 8,
        'lcl_NW_025952766.1_cds_XP_047954044.1_46391__49212': 8,
        'lcl_NC_028637.1_cds_XP_015064474.1_2338__28526': 8,
        'Cc00_g25210.1__49390': 8,
        'rna-gnl_WGS-JAMLDZ_EVM0036354.1__4058': 8,
        'lcl_OX459121.1_cds_CAI9101309.1_14233__43536': 8,
        'lcl_NC_031989.1_cds_XP_019252277.1_2282__49451': 9,
        'Kaladp0011s0492.1.v1.1__63787': 9,
        'lcl_OX459121.1_cds_CAI9101307.1_14231__43536': 9,
        'rna-gnl_WGS-JAMLDZ_EVM0034906.1__4058': 9,
        'lcl_NC_039901.1_cds_XP_027112943.1_15453__13443': 9,
        'rna-gnl_WGS-JAMLDZ_EVM0007978.1__4058': 9,
        'lcl_NC_080155.1_cds_XP_057473373.1_43266__165200': 9,
        'Hma1.2p1_1763F.1_g301360.1__23110': 9,
        'pveT_jg21836.t1__170927': 9,
        'Lj5A64G50.1__105884': 9,
        'lcl_CM028326.1.g5350.t1__4392': 9,
        'rna-gnl_WGS-JAMLDZ_EVM0017794.1__4058': 9,
        'lcl_NC_039905.1_cds_XP_027126787.1_26673__13443': 9,
        'lcl_OX459124.1_cds_CAI9114063.1_26987__43536': 9,
        'lcl_PKPP01007341.1_cds_PWA53652.1_46097__35608': 9,
        'rna-gnl_WGS-JAMLDZ_EVM0014154.1__4058': 9,
        'lcl_OX459124.1_cds_CAI9114065.1_26989__43536': 9,
        'Soltu.DM.12G010750.2__4113': 9,
        'lcl_OU503036.1_cds_CAI9753676.1_1107__56036': 9,
        'lcl_NC_062966.1_cds_XP_047959520.1_10496__49212': 9,
        'lcl_NW_026137468.1_cds_XP_051115997.1_3322__175694': 9,
        'Cc02_g31130.1__49390': 10,
        'Cc00_g28750.1__49390': 10,
        'rna-gnl_WGS-JAMLDZ_EVM0024530.1__4058': 10,
        'Cc03_g08070.1__49390': 10,
        'Anaoc.0007s0766.1__171929': 10,
        'Atru_chr1_3170__47965': 10,
        'FSB015489301__28930': 10,
        'lcl_CM059866.1_cds_KAK0596815.1_1321__4024': 10,
        'FSB011836701__28930': 10,
        'CiLak.01G168500.1__32201': 10,
        'Qurub.10G026200.1__3512': 10,
        'Qurub.10G071600.1__3512': 10,
        'lcl_CM025853.1_cds_KAB1214321.1_15709__262757': 10,
        'CiLak.01G168600.1__32201': 10,
        'lcl_NC_044913.1_cds_XP_030940120.1_44897__97700': 10,
        'Prupe.5G048100.1__3760': 10,
        'Umino00085.1__262084': 10,
        'lcl_NC_024133.1_cds_XP_008243039.1_26638__102107': 10,
        'lcl_NC_084796.1_cds_XP_062119410.1_20878__3486': 10,
        'Umino11360.1__262084': 10,
        'FvH4_3g36530.t1__57918': 10,
        'lcl_NC_083380.1_cds_XP_015881543.2_5733__326968': 10,
        'Kaladp0007s0029.1.v1.1__63787': 10,
        'Blora12794.1__200023': 10,
        'VIT_209s0002g05340.1__29760': 11,
        'lcl_NW_017353139.1_cds_XP_018486826.1_3787__3726': 11,
        'lcl_NC_045129.1_cds_XP_031385398.1_13438__22663': 11,
        'lcl_JAGHRR010000128.1_cds_KAI3418995.1_3930__120290': 11,
        'lcl_NW_026137602.1_cds_XP_051147854.1_28335__175694': 11,
        'Hma1.2p1_0515F.1_g170445.1__23110': 11,
        'lcl_CM028334.1.g444.t1__4392': 11,
        'lcl_OU503038.1_cds_CAI9757760.1_5191__56036': 11,
        'Vadar_g3125.t1__229202': 11,
        'Acc17837.1__3625': 11,
        'pveT_jg21833.t1__170927': 11,
        'Lj5A64T54.1__105884': 11,
        'lcl_NW_026775581.1_cds_XP_059644652.1_44213__4283': 11,
        'Cc03_g10170.1__49390': 11,
        'lcl_CM061492.1_cds_KAK1416820.1_25937__13708': 11,
        'DCAR_016236__4039': 11,
        'lcl_NC_062967.1_cds_XP_047968971.1_19558__49212': 11,
        'TRINITY_DN176594_c2_g1_i2__4102': 11,
        'lcl_NC_045131.1_cds_XP_031399179.1_24312__22663': 11,
        'Eucgr.F02562.1__71139': 11,
        'lcl_NW_006262075.1_cds_XP_006427129.2_25983__85681': 11,
        'lcl_CM059866.1_cds_KAK0598226.1_1319__4024': 11,
        'lcl_NC_052253.1_cds_XP_038693680.1_44458__458696': 11,
        'lcl_NW_011501060.1_cds_XP_011013043.1_5386__75702': 11,
        'lcl_NC_045129.1_cds_XP_031385404.1_13435__22663': 11,
        'Eucgr.F02566.1__71139': 11,
        'Bol038153__3712': 11,
        'Lj1g0016927.1__34305': 11,
        'Ler.1DRT.6g050030.1__41257': 11,
        'Medtr7g090520.1__3880': 11,
        'Vradi08g02800.1__157791': 11,
        'Lj1g0005382.1__34305': 11,
        'Vradi08g02780.1__157791': 11,
        'lcl_CM059866.1_cds_KAK0596259.1_1320__4024': 11,
        'Ptrif.0007s1335.1__37690': 11,
        'Anaoc.0007s0767.1__171929': 11,
        'Gorai.002G057700.1_PACid-26792653__29730': 11,
        'Ptrif.0007s1331.1__37690': 11,
        'lcl_NC_052243.1_cds_XP_038718244.1_22960__458696': 11,
        'lcl_NC_083379.1_cds_XP_015868729.2_4751__326968': 11,
        'lcl_NC_037092.1_cds_XP_024158425.1_33959__74649': 11,
        'Umino20160.1__262084': 11,
        'Umino20144.1__262084': 11,
        'lcl_NC_083601.1_cds_XP_030488467.2_4391__3483': 11,
        'Pparo25676.1__386216': 11,
        'Dglom27811.1__34297': 11,
        'Blora12130.2__200023': 11,
        'lcl_NC_044908.1_cds_XP_030969963.1_22137__97700': 11,
        'lcl_NC_052253.1_cds_XP_038693679.1_44455__458696': 11,
        'Lus10021025_PACid-23182311__4006': 11,
        'lcl_CM058159.1_cds_KAJ9692491.1_17293__103349': 11,
        'Vradi08g02810.1__157791': 11,
        'Ler.1DRT.1g027930.1__41257': 11,
        'Lj5g0008310.1__34305': 11,
        'TRINITY_DN114504_c0_g1_i1__4102': 11,
        'Eucgr.J03023.1__71139': 11,
        'lcl_NC_045131.1_cds_XP_031396688.1_21996__22663': 11,
        'Gorai.011G079100.1_PACid-26810489__29730': 11,
        'Thhalv10011669m__72664': 11,
        'Anaoc.0007s0768.1__171929': 11,
        'lcl_NW_006262075.1_cds_XP_006427132.2_25985__85681': 11,
        'Atru_chr1_3180__47965': 11,
        'Atru_chr1_3175__47965': 11,
        'Potri.001G176000.1__3694': 11,
        'lcl_NC_052244.1_cds_XP_038721828.1_25425__458696': 11,
        'lcl_NC_065574.1_cds_XP_050234356.1_20365__3986': 11,
        'lcl_NC_044908.1_cds_XP_030968259.1_22138__97700': 11,
        'lcl_CM025851.1_cds_KAB1218619.1_7986__262757': 11,
        'CiLak.01G242800.1__32201': 11,
        'lcl_NC_037089.1_cds_XP_024183820.1_11602__74649': 11,
        'lcl_NC_083379.1_cds_XP_048325586.2_4752__326968': 11,
        'lcl_NC_084800.1_cds_XP_062086877.1_35687__3486': 11,
        'Umino20121.1__262084': 11,
        'Dglom27810.1__34297': 11,
        'Tcord01333.1__703396': 11,
        'Blora12127.1__200023': 11,
        'CiLak.02G150500.1__32201': 11,
        'FSB014894201__28930': 11,
        'lcl_NC_044908.1_cds_XP_030968569.1_22136__97700': 11,
        'lcl_PKPP01002761.1_cds_PWA73280.1_26522__35608': 12,
        'rna-gnl_WGS-JAMLDZ_EVM0012381.1__4058': 12,
        'lcl_CM028334.1.g446.t1__4392': 12,
        'Lj5A64T53.1__105884': 12,
        'Hma1.2p1_3098F.1_g332680.1__23110': 12,
        'lcl_NW_026775581.1_cds_XP_059639306.1_44215__4283': 12,
        'Vadar_g3126.t1__229202': 12,
        'Acc12585.1__3625': 12,
        'pveT_jg21835.t1__170927': 12,
        'Pg_S0988.1__4054': 12,
        'DCAR_016235__4039': 12,
        'lcl_OX459124.1_cds_CAI9114062.1_26986__43536': 12,
        'lcl_OU503036.1_cds_CAI9753675.1_1106__56036': 12,
        'lcl_NC_062966.1_cds_XP_047959442.1_10495__49212': 12,
        'lcl_NW_026137602.1_cds_XP_051147853.1_28334__175694': 12,
        'DCAR_007548__4039': 12,
        'lcl_CM061487.1_cds_KAK1434086.1_11004__13708': 12,
        'HanXRQChr09g0275101__4232': 12,
        'Lj5C65T2.1__105884': 12,
        'lcl_NW_026775581.1_cds_XP_059644765.1_44212__4283': 12,
        'lcl_NC_080143.1_cds_XP_057503960.1_20127__165200': 12,
        'pveT_jg21832.t1__170927': 12,
        'lcl_CM028334.1.g443.t1__4392': 12,
        'rna-gnl_WGS-JAMLDZ_EVM0018858.1__4058': 12,
        'SMEL_006g260960.1.01__4111': 12,
        'lcl_NC_039905.1_cds_XP_027126923.1_26675__13443': 12,
        'lcl_OU503036.1_cds_CAI9753673.1_1104__56036': 12,
        'MSTRG.22708.1__82927': 12,
        'lcl_NC_062966.1_cds_XP_047960078.1_10492__49212': 12,
        'lcl_NW_026137587.1_cds_XP_051135949.1_16701__175694': 12,
        'lcl_NW_026137587.1_cds_XP_051135915.1_16702__175694': 12,
        'lcl_CM008454.1_cds_PHT32940.1_29277__33114': 13,
        'lcl_CM008444.1_cds_PHT55775.1_4261__33114': 13,
        'TRINITY_DN232929_c4_g6_i3__4102': 13,
        'lcl_CM008445.1_cds_PHT51893.1_6355__33114': 13,
        'lcl_NC_083337.1_cds_XP_060201071.1_373__112863': 13,
        'lcl_NC_083338.1_cds_XP_060200949.1_5353__112863': 13,
        'Soltu.DM.01G002190.1__4113': 13,
        'SMEL_001g152420.1.01__4111': 13,
        'lcl_NC_028638.1_cds_XP_015067085.1_5485__28526': 13,
        'TRINITY_DN725607_c0_g1_i1__189803': 13,
        'Vadar_g43141.t1__229202': 13,
        'KT390173__DPS__etoposide_biosynthesis__93608': 13,
        'lcl_NW_026775571.1_cds_XP_059664661.1_3040__4283': 13,
        'AF417859__AOP3__glucosinolate_biosynthesis__3702': 13,
        'AF417858__AOP2__glucosinolate_biosynthesis__3702': 13,
        'lcl_CM042048.1_cds_KAI3759661.1_7640__4217': 13,
        'lcl_NW_026137594.1_cds_XP_051141059.1_21911__175694': 13,
        'lcl_NC_083379.1_cds_XP_048320867.2_1719__326968': 13,
        'Umino20164.1__262084': 13,
        'Cc08_g04750.1__49390': 13,
        'Lus10023024_PACid-23146145__4006': 13,
        'lcl_NC_045131.1_cds_XP_031396066.1_23250__22663': 13,
        'lcl_NC_052241.1_cds_XP_038713292.1_19041__458696': 13,
        'Gobar.A11G255200.1__3634': 13,
        'Medtr4g011690.1__3880': 13,
        'Ca_13466__3827': 13,
        'Ler.1DRT.4g003660.1__41257': 13,
        'Ca_13467__3827': 13,
        'Lj3g0013768.1__34305': 13,
        'lcl_OX451740.1_cds_CAI8613451.1_28042__3906': 13,
        'Lj1g0004205.1__34305': 13,
        'Vradi0154s00060.1__157791': 13,
        'lcl_NC_052245.1_cds_XP_038722105.1_29001__458696': 13,
        'lcl_NC_065574.1_cds_XP_050231702.1_19581__3986': 13,
        'Lus10005516_PACid-23158832__4006': 13,
        'lcl_NC_052393.1_cds_XP_008784635.1_5049__42345': 13,
        'Vadar_g3127.t1__229202': 13,
        'VIT_209s0002g05280.1__29760': 13,
        'Acc17838.1__3625': 13,
        'lcl_NW_026775581.1_cds_XP_059643089.1_44218__4283': 13,
        'lcl_CM028334.1.g447.t1__4392': 13,
        'Lj9C505T3.1__105884': 13,
        'Pg_S0213.71__4054': 13,
        'DCAR_016234__4039': 13,
        'lcl_CM027379.1_cds_KAF9612984.1_11656__261450': 13,
        'lcl_PKPP01007341.1_cds_PWA53651.1_46096__35608': 13,
        'Hma1.2p1_0515F.1_g170455.1__23110': 13,
        'Ah1G51780.1__81970': 13,
        'Cvari18296.1__869952': 13,
        'Umino20166.1__262084': 13,
        'lcl_NC_083379.1_cds_XP_015868730.2_4749__326968': 13,
        'lcl_NC_044908.1_cds_XP_030973050.1_22133__97700': 13,
        'CiLak.01G242700.1__32201': 13,
        'lcl_CM025851.1_cds_KAB1218616.1_7983__262757': 13,
        'Prupe.3G075700.1__3760': 13,
        'Dglom27813.1__34297': 13,
        'Blora09835.1__200023': 13,
        'lcl_NW_006262075.1_cds_XP_006427127.2_25981__85681': 13,
        'Potri.001G176200.1__3694': 13,
        'lcl_NC_065574.1_cds_XP_050231462.1_20561__3986': 13,
        'lcl_JAGHRR010000039.1_cds_KAI3437248.1_15716__120290': 13,
        'Anaoc.0007s0764.1__171929': 13,
        'lcl_CM059866.1_cds_KAK0595174.1_1315__4024': 13,
        'OMO54201__210143': 13,
        'lcl_NC_052253.1_cds_XP_038694945.1_44457__458696': 13,
        'Lus10016659_PACid-23143917__4006': 13,
        'pveT_jg11526.t1__170927': 13,
        'lcl_NC_083601.1_cds_XP_060957775.1_4387__3483': 13,
        'lcl_OX459123.1_cds_CAI9108625.1_21549__43536': 13,
        'Kaladp0008s0237.1.v1.1__63787': 13,
        'Soltu.DM.06G023510.1__4113': 13,
        'rna-gnl_WGS-JAMLDZ_EVM0014967.1__4058': 13,
        'lcl_OU503036.1_cds_CAI9753677.1_1108__56036': 13,
        'MSTRG.11101.1__82927': 13,
        'lcl_NC_062966.1_cds_XP_047958340.1_10499__49212': 13,
        'lcl_NC_062965.1_cds_XP_047956864.1_861__49212': 13,
        'Mba05_g08120.1.v1.1__52838': 13,
        'CALSI_Maker00003989__746888': 13,
        'Mba03_g21520.1.v1.1__52838': 13,
        'evm.model.AsparagusV1_08.662__4686': 13,
        'Ma04_t25080.1__4641': 13,
        'lcl_CP136891.1_cds_WOK98631.1_7322__4628': 13,
        'Mba05_g08100.1.v1.1__52838': 13,
        'Mtr04_g32474.1__320322': 13,
        'lcl_CP136893.1_cds_WOL03332.1_12022__4628': 13,
        'CALSI_Maker00009131__746888': 13,
        'ORUFI05G28750.1__4529': 13,
        'ORUFI11G09650.1__4529': 13,
        'lcl_NC_052396.1_cds_XP_008794218.1_12468__42345': 13,
        'Spipo4G0051700__29656': 13,
        'lcl_CP136890.1_cds_WOK92059.1_750__4628': 13,
        'CALSI_Maker00050580__746888': 13,
        'evm.model.AsparagusV1_01.591__4686': 13,
        'ASH_rna4294__1088818': 13,
        'strangu_020525-RA__38733': 13,
        'lcl_NC_053024.1_cds_XP_048567703.1_14633__4572': 14
    }

    result = seq_id_to_landscape_idx(resolved_landscape)
    assert result == expected 


def test_landscape_meta_info():
    manual_candidates = {
        "lcl_NC_084852.1_cds_XP_061960601.1_3168__3691" : "2ODD37",
        'VIT_204s0008g04920.2__29760' : '2ODD37',
        'Ptrif.0001s0318.1__37690' : '2ODD37',
        'Atru_chr7_2342__47965' : '2ODD37',
        'lcl_NW_019168159.1_cds_XP_022728803.1_47147__66656' : '2ODD37',
        'lcl_NC_045134.1_cds_XP_031407529.1_33286__22663' : '2ODD37',
        'GWHTACBH000860__413952' : '2ODD37',
        'Potri.006G248000.1__3694' : '2ODD37',
        'lcl_NC_065570.1_cds_XP_050209042.1_2787__3986' : '2ODD37',
        'Lus10008097_PACid-23169790__4006' : '2ODD37', 
        'Acc09929.1__3625': "minor_2ODD_cluster",
        'lcl_NC_031989.1_cds_XP_019252277.1_2282__49451': "2ODD38",
        'Kaladp0011s0492.1.v1.1__63787': "2ODD38",
        'lcl_OX459121.1_cds_CAI9101307.1_14231__43536': "2ODD38",
        'rna-gnl_WGS-JAMLDZ_EVM0034906.1__4058': "2ODD38",
        'lcl_NC_039901.1_cds_XP_027112943.1_15453__13443': "2ODD38",
        'rna-gnl_WGS-JAMLDZ_EVM0007978.1__4058': "2ODD38", 
        'Bol038153__3712': "2ODD39",
        "lcl_NC_053024.1_cds_XP_048567703.1_14633__4572": "unresolved"
        }

    for leaf in t:
        if leaf.name in manual_candidates:
            leaf.props["two_odd_id"] = "candidate"
    _, ingroup_headers = split_seqs_by_2ODD_membership(seq_to_2ODD_id=seq_to_2ODD_id, tree=t)

    ingroup_headers = ingroup_headers - set(manual_candidates.keys())
    landscape_unresolved = get_landscape(t)
    resolved_landscape = resolve_candidates_in_landscape(landscape=landscape_unresolved, 
                                                            dist_dict=dist_dict, 
                                                            threshold=2)
    
    expected = {
        0: {'two_odd_id': '2ODD37',
        'percentage_of_ingroup_2ODD': 1.0,
        'num_sequences': 48,
        'plant_groups': ['Dicots'],
        'contains_candidates': True},
        1: {'two_odd_id': 'minor_2ODD_cluster',
        'percentage_of_ingroup_2ODD': None,
        'num_sequences': 4,
        'plant_groups': ['Dicots'],
        'contains_candidates': True},
        2: {'two_odd_id': '2ODD41',
        'percentage_of_ingroup_2ODD': 0.078,
        'num_sequences': 8,
        'plant_groups': ['Dicots'],
        'contains_candidates': False},
        3: {'two_odd_id': 'minor_2ODD_cluster',
        'percentage_of_ingroup_2ODD': None,
        'num_sequences': 1,
        'plant_groups': ['Dicots'],
        'contains_candidates': False},
        4: {'two_odd_id': '2ODD40',
        'percentage_of_ingroup_2ODD': 0.667,
        'num_sequences': 20,
        'plant_groups': ['Dicots'],
        'contains_candidates': False},
        5: {'two_odd_id': 'minor_2ODD_cluster',
        'percentage_of_ingroup_2ODD': None,
        'num_sequences': 9,
        'plant_groups': ['Dicots'],
        'contains_candidates': False},
        6: {'two_odd_id': '2ODD40',
        'percentage_of_ingroup_2ODD': 0.033,
        'num_sequences': 1,
        'plant_groups': ['Dicots'],
        'contains_candidates': False},
        7: {'two_odd_id': 'minor_2ODD_cluster',
        'percentage_of_ingroup_2ODD': None,
        'num_sequences': 3,
        'plant_groups': ['Dicots'],
        'contains_candidates': False},
        8: {'two_odd_id': '2ODD40',
        'percentage_of_ingroup_2ODD': 0.3,
        'num_sequences': 9,
        'plant_groups': ['Dicots'],
        'contains_candidates': False},
        9: {'two_odd_id': '2ODD38',
        'percentage_of_ingroup_2ODD': 0.326,
        'num_sequences': 21,
        'plant_groups': ['Dicots'],
        'contains_candidates': True},
        10: {'two_odd_id': 'minor_2ODD_cluster',
        'percentage_of_ingroup_2ODD': None,
        'num_sequences': 24,
        'plant_groups': ['Dicots'],
        'contains_candidates': False},
        11: {'two_odd_id': '2ODD39',
        'percentage_of_ingroup_2ODD': 1.0,
        'num_sequences': 79,
        'plant_groups': ['Dicots'],
        'contains_candidates': True},
        12: {'two_odd_id': '2ODD38',
        'percentage_of_ingroup_2ODD': 0.674,
        'num_sequences': 31,
        'plant_groups': ['Dicots'],
        'contains_candidates': False},
        13: {'two_odd_id': '2ODD41',
        'percentage_of_ingroup_2ODD': 0.922,
        'num_sequences': 95,
        'plant_groups': ['Dicots', 'Monocots'],
        'contains_candidates': False},
        14: {'two_odd_id': 'unresolved',
        'percentage_of_ingroup_2ODD': None,
        'num_sequences': 1,
        'plant_groups': ['Monocots'],
        'contains_candidates': True}
        }

    result = landscape_meta_info(resolved_landscape, major_minor_2ODDs_dict,ingroup_headers, manual_candidates)
    assert result == expected

