from two_odd_annotator.utils.logging import LOG_FILENAME, init_log, log_line


def test_init_log_creates_file_and_header(tmp_path, monkeypatch):
    # Make timestamps deterministic for the header check
    fake_time = "2026-03-01T12:00:00+00:00"
    monkeypatch.setattr(
        "two_odd_annotator.utils.logging._now_iso",
        lambda: fake_time,
    )

    log_path = init_log(str(tmp_path))

    assert log_path.name == LOG_FILENAME
    assert log_path.is_file()

    content = log_path.read_text().splitlines()
    assert content[0] == f"=== Run started {fake_time} ==="


def test_log_line_appends_message(tmp_path, monkeypatch):
    # Initialise log file first
    log_path = init_log(str(tmp_path))

    fake_time = "2026-03-01T12:34:56+00:00"
    monkeypatch.setattr(
        "two_odd_annotator.utils.logging._now_iso",
        lambda: fake_time,
    )

    log_line(logfile_path=log_path, message="Hello world")

    lines = (tmp_path / LOG_FILENAME).read_text().splitlines()
    # Last line should be the message we just logged
    assert lines[-1] == f"{fake_time} - Hello world" 
