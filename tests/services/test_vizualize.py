import json

import pandas as pd


def test_load_presence_major_only(tmp_path):
    from two_odd_annotator.services.analyze import _load_presence_by_taxid

    major_minor = {
        "major_2ODDs": {
            "2ODD01": ["ING_A__3702"],
            "2ODD02": ["ING_B__4530"],
        },
        "minor_2ODDs": {
            "2ODD_minor14": ["MIN_A__3702"],
        },
    }
    major_minor_path = tmp_path / "major_minor.json"
    major_minor_path.write_text(json.dumps(major_minor))

    # annotation TSV includes mixed assignments
    df = pd.DataFrame(
        {
            "candidate": [
                "C1__3702",
                "C2__3702",
                "C3__4530",
                "C4__4530",
            ],
            "cluster_two_odd_id": [
                "2ODD01",  # keep
                "2ODD_minor14",  # drop (minor)
                "candidates_only",  # drop (candidate-only cluster)
                "minor_2ODD_cluster",  # drop (placeholder)
            ],
        }
    )
    annotation_tsv = tmp_path / "annotation_results.tsv"
    df.to_csv(annotation_tsv, sep="\t", index=False)

    presence, major_ids = _load_presence_by_taxid(
        output_dir=tmp_path,
        major_minor_json=major_minor_path,
        annotation_tsv=annotation_tsv,
    )

    assert major_ids == {"2ODD01", "2ODD02"}
    assert presence[3702] == {"2ODD01"}
    assert presence[4530] == {"2ODD02"}


def test_make_ultrametric_biophylo_equalizes_tip_depths():
    from io import StringIO

    from Bio import Phylo

    from two_odd_annotator.services.analyze import _make_ultrametric_biophylo

    # Non-ultrametric tree: tips have different depths.
    t = Phylo.read(StringIO("(A:1,(B:1,C:2):1);"), "newick")
    depths_before = t.depths()
    tip_depths_before = [depths_before[x] for x in t.get_terminals()]
    assert len(set(tip_depths_before)) > 1

    _make_ultrametric_biophylo(t)
    depths_after = t.depths()
    tip_depths_after = [round(depths_after[x], 6) for x in t.get_terminals()]
    assert len(set(tip_depths_after)) == 1


def test_build_display_labels_appends_species_counts():
    from two_odd_annotator.services.analyze import _build_display_labels

    labels = {1: "OrderA", 2: "OrderB"}
    counts = {1: 12, 2: 1}
    out = _build_display_labels(labels, counts)

    assert out[1] == "OrderA (12)"
    assert out[2] == "OrderB (1)"


def test_classify_plant_group_examples():
    from two_odd_annotator.services.analyze import _classify_plant_group

    assert _classify_plant_group(["cellular organisms", "viridiplantae", "liliopsida"]) == "Monocots"
    assert _classify_plant_group(["viridiplantae", "eudicotyledons"]) == "Dicots"
    assert _classify_plant_group(["viridiplantae", "polypodiopsida"]) == "Ferns"
    assert _classify_plant_group(["viridiplantae", "acrogymnospermae"]) == "Gymnosperms"


def test_sort_biophylo_tree_by_terminal_order_reorders_children():
    from io import StringIO

    from Bio import Phylo

    from two_odd_annotator.services.analyze import _sort_biophylo_tree_by_terminal_order

    # Start with a tree whose terminals are in A,B,C order.
    t = Phylo.read(StringIO("(1:1,(2:1,3:1):1);"), "newick")
    before = [x.name for x in t.get_terminals()]
    assert before == ["1", "2", "3"]

    # Force desired order 3,2,1.
    _sort_biophylo_tree_by_terminal_order(t, [3, 2, 1])
    after = [x.name for x in t.get_terminals()]
    assert after == ["3", "2", "1"]


def test_species_taxids_for_analysis_drops_char_only_species(tmp_path):
    from two_odd_annotator.services.analyze import _species_taxids_for_analysis

    # taxid 111: characterized-only => should be dropped
    # taxid 222: characterized + at least one non-characterized => keep
    # taxid 333: non-characterized only => keep
    fasta = tmp_path / "annotation.fasta"
    fasta.write_text(
        ">CHAR__X__Y__111\nAAAA\n"
        ">CHAR2__X__Y__111\nAAAA\n"
        ">CHAR__X__Y__222\nAAAA\n"
        ">CAND__222\nAAAA\n"
        ">CAND__333\nAAAA\n"
    )

    kept, dropped = _species_taxids_for_analysis(fasta)
    assert 111 in dropped
    assert 111 not in kept
    assert 222 in kept
    assert 333 in kept


def test_two_odd_xtick_labels_appends_associated_functions():
    from two_odd_annotator.services.analyze import _two_odd_xtick_labels

    major_char = {
        "2ODD01": {"associated_functions": ["D4H"]},
        "2ODD02": {"associated_functions": ["E8", "GAME40"]},
        "2ODD03": {"associated_functions": []},
    }

    out = _two_odd_xtick_labels(["2ODD01", "2ODD02", "2ODD03", "2ODDXX"], major_char)
    assert out[0] == "(D4H) 2ODD01"
    assert out[1] == "(E8,GAME40) 2ODD02"
    assert out[2] == "2ODD03"
    assert out[3] == "2ODDXX"
