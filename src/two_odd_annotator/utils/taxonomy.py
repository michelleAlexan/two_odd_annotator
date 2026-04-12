from __future__ import annotations

from typing import Dict

from ete4.ncbi_taxonomy import NCBITaxa


_ncbi = NCBITaxa()

# Handle a few known spelling/spacing variants that do not resolve
# cleanly via the NCBI name translator.
_ALIASES: dict[str, str] = {
    "Echinochloa crus galli": "Echinochloa crus-galli",
}


def map_scientific_notation_to_tax_id(
    scientific_name: str,
    raise_on_error: bool = True,
) -> Dict[str, int]:
    """Map a scientific species name to its NCBI taxonomic ID.

    Parameters
    ----------
    scientific_name : str
        Species name in "Genus species" form. Underscores are treated as
        spaces for convenience.
    raise_on_error : bool, default True
        If True, raise ValueError when no taxid can be found. If False,
        return an empty dict.

    Returns
    -------
    dict[str, int]
        A one-element mapping ``{resolved_name: taxid}``.
    """

    name = scientific_name.replace("_", " ").strip()
    name = _ALIASES.get(name, name)

    try:
        name_to_ids = _ncbi.get_name_translator([name])
    except Exception as exc:  # pragma: no cover - defensive
        if raise_on_error:
            raise ValueError(
                f"Failed to query NCBI taxonomy for: {scientific_name}"
            ) from exc
        return {}

    if not name_to_ids:
        if raise_on_error:
            raise ValueError(
                f"Could not map species name to tax ID: {scientific_name}"
            )
        return {}

    resolved_name, ids = next(iter(name_to_ids.items()))
    taxid = int(ids[0])
    return {resolved_name: taxid}
