from __future__ import annotations

from dataclasses import asdict, is_dataclass
from datetime import datetime, timezone
from email import message
from pathlib import Path
from typing import Any, Mapping


LOG_FILENAME = "run.log"


def _now_iso() -> str:
    """Return current UTC time as an ISO 8601 string."""

    return datetime.now(timezone.utc).isoformat()


def init_log(output_dir: str, logfile_name: str = LOG_FILENAME) -> Path:
    """Create (or append to) a simple text log file.

    - Ensures the output directory exists.
    - Writes a header with the current datetime.
    - Optionally dumps the CLI arguments / run configuration.
    """

    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    log_path = out_dir / logfile_name
    if log_path.is_file():
        # If the log already exists, it is overwritten to start fresh for each run, 
        # but you could change this behavior to keep a cumulative log.
        log_path.unlink()

    with log_path.open("w") as f:
        f.write(f"=== Run started { _now_iso() } ===\n")
        f.write("\n")
        print(f"Initialized log at {log_path}")
    return log_path


def log_line(logfile_path: str, message: str) -> None:
    """Append a timestamped message to the run log in output_dir."""

    log_path = Path(logfile_path)
    with log_path.open("a") as f:
        f.write(f"{_now_iso()} - {message}\n")
    print(message)  # Also print to console for real-time feedback

__all__ = ["LOG_FILENAME", "init_log", "log_line"]



