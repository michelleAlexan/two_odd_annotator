import json


def test_state_applies_species_name_mapping(monkeypatch, tmp_path):
    # Create an input FASTA that would otherwise fail mapping.
    input_dir = tmp_path / "in"
    input_dir.mkdir()
    fasta = input_dir / "Gynostemma_pentaphyllum1.pep.fa"
    fasta.write_text(">seq1\nACGT\n")

    # Provide a mapping that removes the trailing "1".
    mapping = {"Gynostemma pentaphyllum1": "Gynostemma pentaphyllum"}
    mapping_path = tmp_path / "sp_map.json"
    mapping_path.write_text(json.dumps(mapping))

    # Monkeypatch taxonomy lookup to validate the mapped name is used.
    import two_odd_annotator.pipeline.state as state_mod

    seen = {}

    def fake_map(name: str, raise_on_error: bool = True):
        seen["name"] = name
        return {name: 999}

    monkeypatch.setattr(state_mod, "map_scientific_notation_to_tax_id", fake_map)

    from two_odd_annotator.pipeline.state import State

    out_dir = tmp_path / "out"
    s = State(input_path=input_dir, output_base_dir=out_dir, sp_name_mapping_path=mapping_path)

    assert seen["name"] == "Gynostemma pentaphyllum"
    # Ensure it still initialised a species subdir.
    assert len(s.metadata) == 1
