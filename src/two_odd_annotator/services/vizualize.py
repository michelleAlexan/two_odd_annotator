from pathlib import Path
from typing import Any


def run(output_dir: Path, config: dict[str, Any]) -> None:
    """Placeholder visualization step.

    This function exists so that the ``visualize`` pipeline step can be
    invoked from the CLI without raising an AttributeError. The actual
    plotting/visualization logic can be implemented here in the future.
    """

    print(f"Visualization step is not implemented yet. Output dir: {output_dir}")
