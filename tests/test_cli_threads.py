import sys


def test_cli_threads_forwarded_to_runner(monkeypatch, tmp_path):
    import two_odd_annotator.cli as cli

    captured = {}

    class DummyRunner:
        def __init__(self, *args, **kwargs):
            captured.update(kwargs)

        def run(self):
            return None

    monkeypatch.setattr(cli, "Runner", DummyRunner)

    out_dir = tmp_path / "out"
    argv = [
        "annodd",
        "-i",
        str(tmp_path),
        "-o",
        str(out_dir),
        "--threads",
        "7",
    ]
    monkeypatch.setattr(sys, "argv", argv)

    cli.main()

    assert captured.get("threads") == 7


def test_runner_threads_override_config(tmp_path):
    from pathlib import Path

    from two_odd_annotator.pipeline.runner import Runner

    test_data_path = Path(__file__).parents[0] / "data"
    test_config_path = Path(__file__).parents[0] / "config" / "test_config.yml"

    output_dir = tmp_path / "results"

    runner = Runner(
        input_path=str(test_data_path),
        output_base_dir=str(output_dir),
        config_path=str(test_config_path),
        threads=5,
        step="filter_seq_sim",  
    )

    assert runner.config["parameters"]["threads"] == 5
