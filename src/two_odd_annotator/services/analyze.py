from __future__ import annotations

from pathlib import Path
from typing import Any
import json
import re

import pandas as pd
from Bio import SeqIO
from ete4 import NCBITaxa

from two_odd_annotator.constants import (
    ANNOTATION_FASTA,
    ANNOTATION_TSV,
)
from two_odd_annotator.utils.logging import LOG_FILENAME, log_line


GROUP_COLORS: dict[str, str] = {
    "Algae": "#574104",
    "Lycophytes": "#ab730c",
    "Liverworts": "#faaf00",
    "Ferns": "#c1d717",
    "Mosses": "#89be86",
    "Gymnosperms": "#076247",
    "Basal Angiosperms": "#3BA0BC",
    "Monocots": "#7502d9",
    "Dicots": "#edc5ec",
}


def _classify_plant_group(named_lineage_lower: list[str]) -> str | None:
    """Classify a lineage into a high-level plant group.

    Expects lower-cased scientific names from an NCBI named lineage.
    Returns None when no mapping applies.
    """

    lineage = named_lineage_lower

    if "zygnematophyceae" in lineage:
        return "Algae"
    if "lycopodiopsida" in lineage:
        return "Lycophytes"
    if "polypodiopsida" in lineage:
        return "Ferns"
    if "marchantiophyta" in lineage:
        return "Liverworts"
    if ("bryophyta" in lineage) or ("anthocerotophyta" in lineage):
        return "Mosses"
    if any(
        x in lineage
        for x in [
            "amborellales",
            "nymphaeales",
            "austrobaileyales",
            "magnoliidae",
        ]
    ):
        return "Basal Angiosperms"
    if "acrogymnospermae" in lineage:
        return "Gymnosperms"
    if "liliopsida" in lineage:
        return "Monocots"
    if any(x in lineage for x in ["eudicotyledons", "magnoliopsida", "mesangiospermae"]):
        return "Dicots"

    return None


def _plant_group_for_taxid(ncbi: NCBITaxa, taxid: int) -> str | None:
    """Return the plant group name for a taxid, or None if not classifiable."""

    try:
        lineage = ncbi.get_lineage(int(taxid))
        names = ncbi.get_taxid_translator(lineage)
        named_lineage_lower = [str(names.get(t, "")).lower() for t in lineage]
    except Exception:
        return None

    return _classify_plant_group(named_lineage_lower)


def _label_color_for_taxid(ncbi: NCBITaxa, taxid: int) -> str | None:
    """Return a hex color for a taxid label, or None for default styling."""

    group = _plant_group_for_taxid(ncbi, taxid)
    if not group:
        return None
    return GROUP_COLORS.get(group)


def _two_odd_xtick_labels(
    two_odd_ids: list[str],
    major_char_info: dict[str, Any] | None,
) -> list[str]:
    """Build display labels for 2ODD IDs with associated functions."""

    if not major_char_info:
        return list(two_odd_ids)

    labels: list[str] = []
    for tid in two_odd_ids:
        info = major_char_info.get(tid, {}) if isinstance(major_char_info, dict) else {}
        funcs = info.get("associated_functions") if isinstance(info, dict) else None
        if isinstance(funcs, list):
            funcs = [str(x) for x in funcs if x is not None and str(x).strip()]
        else:
            funcs = []

        if funcs:
            labels.append(f"({','.join(funcs)}) {tid}")
        else:
            labels.append(tid)
    return labels


def _build_display_labels(
    row_labels: dict[int, str],
    species_counts: dict[int, int] | None,
) -> dict[int, str]:
    """Return taxid->label strings, optionally appending " (N)" for species counts."""

    if not species_counts:
        return dict(row_labels)

    out: dict[int, str] = {}
    for taxid, name in row_labels.items():
        n = int(species_counts.get(taxid, 0))
        if n > 0:
            out[taxid] = f"{name} ({n})"
        else:
            out[taxid] = name
    return out


def _extract_taxid_from_seq_id(seq_id: str) -> int | None:
    """Extract NCBI taxid from a sequence ID.

    Expected convention in this project: sequence IDs end with "__<taxid>".
    Returns None if no trailing integer taxid is found.
    """

    if "__" in seq_id:
        last = seq_id.split("__")[-1]
        if last.isdigit():
            return int(last)
    m = re.search(r"(\d+)$", seq_id)
    if m:
        return int(m.group(1))
    return None


def _is_char_bait_sequence(seq_id: str) -> bool:
    """Return True if seq_id matches the characterized-bait ID pattern.

    This mirrors the annotation logic in services/annotate.py:
    characterized bait sequences have exactly 4 "__"-separated fields.
    """

    return len(seq_id.split("__")) == 4


def _species_taxids_for_analysis(annotation_fasta: Path) -> tuple[list[int], list[int]]:
    """Return (kept_taxids, dropped_taxids) based on dataset availability.

    For each taxid that appears in the annotation FASTA as a characterized bait,
    require at least one *non-characterized* sequence with the same taxid.
    If not present, drop that taxid from downstream analysis.

    Rationale: some characterized sequences can be included even when the full
    species dataset was not part of the run, which would bias presence/absence.
    """

    char_counts: dict[int, int] = {}
    non_char_counts: dict[int, int] = {}

    for rec in SeqIO.parse(str(annotation_fasta), "fasta"):
        tid = _extract_taxid_from_seq_id(rec.id)
        if tid is None:
            continue
        if _is_char_bait_sequence(rec.id):
            char_counts[tid] = int(char_counts.get(tid, 0)) + 1
        else:
            non_char_counts[tid] = int(non_char_counts.get(tid, 0)) + 1

    all_taxids = sorted(set(char_counts) | set(non_char_counts))
    dropped_set = {tid for tid in all_taxids if char_counts.get(tid, 0) > 0 and non_char_counts.get(tid, 0) == 0}
    dropped = sorted(dropped_set)
    kept = [tid for tid in all_taxids if tid not in dropped_set]
    return kept, dropped


def _load_presence_by_taxid(
    output_dir: Path,
    major_minor_json: Path,
    annotation_tsv: Path,
) -> tuple[dict[int, set[str]], set[str]]:
    """Return {taxid -> set(major 2ODD IDs present)} and the major 2ODD ID set.

    Heatmaps intentionally exclude:
    - minor 2ODD IDs (from the JSON)
    - candidate-only clusters (e.g. "candidates_only")
    - minor cluster placeholder IDs (e.g. "minor_2ODD_cluster")

    Candidate sequences are still counted if they were assigned to a *major* 2ODD ID.
    """

    presence: dict[int, set[str]] = {}
    major_two_odd_ids: set[str] = set()

    major_minor = json.load(open(major_minor_json))
    for two_odd_id, seq_ids in major_minor.get("major_2ODDs", {}).items():
        major_two_odd_ids.add(two_odd_id)
        for seq_id in seq_ids:
            tid = _extract_taxid_from_seq_id(seq_id)
            if tid is None:
                continue
            presence.setdefault(tid, set()).add(two_odd_id)

    # candidates: cluster-level 2ODD assignment from annotation_results.tsv
    anno_df = pd.read_csv(annotation_tsv, sep="\t")
    if "candidate" not in anno_df.columns:
        raise ValueError(f"Expected column 'candidate' in {annotation_tsv.name}")
    if "cluster_two_odd_id" in anno_df.columns:
        col = "cluster_two_odd_id"
    elif "two_odd_id" in anno_df.columns:
        col = "two_odd_id"
    else:
        raise ValueError(
            "Expected a 2ODD assignment column ('cluster_two_odd_id' or 'two_odd_id') in annotation results"
        )

    for _, row in anno_df.iterrows():
        cand = str(row["candidate"])
        tid = _extract_taxid_from_seq_id(cand)
        if tid is None:
            continue
        two_odd_id = row.get(col)
        if two_odd_id is None or (isinstance(two_odd_id, float) and pd.isna(two_odd_id)):
            continue
        two_odd_id = str(two_odd_id)
        if two_odd_id in major_two_odd_ids:
            presence.setdefault(tid, set()).add(two_odd_id)

    return presence, major_two_odd_ids


def _rank_taxid(ncbi: NCBITaxa, taxid: int, rank: str) -> int | None:
    if rank == "species":
        return taxid

    lineage = ncbi.get_lineage(taxid)
    ranks = ncbi.get_rank(lineage)
    for tid in lineage:
        if ranks.get(tid) == rank:
            return tid
    return None


def _topology_order(ncbi: NCBITaxa, taxids: list[int]) -> list[int]:
    """Return taxids ordered by an NCBI taxonomy topology."""
    if not taxids:
        return []
    if len(taxids) == 1:
        return taxids
    t = ncbi.get_topology(taxids, annotate=False)
    return [int(leaf.name) for leaf in t.leaves() if str(leaf.name).isdigit()]


def _make_ultrametric_biophylo(tree) -> None:
    """Transform a Bio.Phylo tree into an ultrametric tree in-place.

    This adjusts terminal branch lengths so that all root→tip distances are equal.
    Missing branch lengths are treated as 1.0.
    """

    # Ensure branch lengths exist (Bio.Phylo uses None by default)
    for clade in tree.find_clades(order="preorder"):
        if clade.branch_length is None:
            clade.branch_length = 1.0

    depths = tree.depths()  # distances from root to each clade
    terminals = tree.get_terminals()
    if not terminals:
        return

    max_depth = max(depths.get(t, 0.0) for t in terminals)
    for t in terminals:
        d = depths.get(t, 0.0)
        # extend terminal branch so that all tips end at max_depth
        t.branch_length = float(t.branch_length) + float(max_depth - d)


def _sort_biophylo_tree_by_terminal_order(tree, desired_taxid_order: list[int]) -> None:
    """Sort Bio.Phylo tree clade order to match a desired terminal ordering.

    This reorders children clades in-place so the leaf (tip) order in the drawn
    tree is consistent with the matrix row order.
    """

    if not desired_taxid_order:
        return

    order = {str(t): i for i, t in enumerate(desired_taxid_order)}
    default_idx = len(order) + 1

    def _clade_min_idx(clade) -> int:
        mins: list[int] = []
        for term in clade.find_clades(terminal=True):
            if term.name is None:
                continue
            mins.append(order.get(str(term.name), default_idx))
        return min(mins) if mins else default_idx

    def _sort_clade(clade) -> None:
        if not getattr(clade, "clades", None):
            return
        for child in clade.clades:
            _sort_clade(child)
        clade.clades.sort(key=_clade_min_idx)

    _sort_clade(tree.root)


def run(output_dir: Path, config: dict[str, Any]) -> None:
    log_path = str(Path(output_dir) / LOG_FILENAME)

    output_dir = Path(output_dir)
    annotation_tsv = output_dir / ANNOTATION_TSV
    annotation_fasta = output_dir / ANNOTATION_FASTA

    if not annotation_tsv.exists():
        log_line(log_path, f"[ANALYZE] Skipping: {ANNOTATION_TSV} not found in {output_dir}")
        return
    if not annotation_fasta.exists():
        log_line(log_path, f"[ANALYZE] Skipping: {ANNOTATION_FASTA} not found in {output_dir}")
        return

    # Prefer new config section name, but keep backwards compatibility.
    analyze_cfg = config.get("analyze")
    if analyze_cfg is None:
        analyze_cfg = config.get("visualize", {})
    analyze_cfg = analyze_cfg or {}

    rank = str(analyze_cfg.get("heatmap_taxa_rank", "order")).lower()
    if rank not in {"species", "genus", "family", "order"}:
        raise ValueError("analyze.heatmap_taxa_rank must be one of: species, genus, family, order")

    major_minor_json = Path(config["annotate"]["major_minor_2ODD_ids"])
    major_char_info_path = Path(config["annotate"]["major_2ODDs_functional_characterization"])

    major_char_info: dict[str, Any] | None = None
    if major_char_info_path.exists():
        try:
            major_char_info = json.load(open(major_char_info_path))
        except Exception as e:
            log_line(log_path, f"[ANALYZE] Could not read major characterization JSON ({major_char_info_path}): {e}")
            major_char_info = None
    else:
        log_line(
            log_path,
            f"[ANALYZE] Major characterization JSON not found at {major_char_info_path}; x labels will be IDs only",
        )

    log_line(log_path, f"[ANALYZE] Building presence heatmap at rank={rank}")

    # gather species taxids from annotation fasta and filter out taxa that only
    # appear as characterized baits without any dataset sequences.
    species_taxids, dropped_species_taxids = _species_taxids_for_analysis(annotation_fasta)
    if not species_taxids:
        log_line(log_path, "[ANALYZE] No taxids found in annotation.fasta headers; skipping")
        return
    if dropped_species_taxids:
        log_line(
            log_path,
            f"[ANALYZE] Dropping {len(dropped_species_taxids)} species that appear only as characterized baits (no dataset sequences)",
        )

    presence_by_tid, two_odd_ids = _load_presence_by_taxid(
        output_dir=output_dir,
        major_minor_json=major_minor_json,
        annotation_tsv=annotation_tsv,
    )
    two_odd_ids_sorted = sorted(two_odd_ids)

    # species-level matrix (base for aggregation)
    base = pd.DataFrame(
        0,
        index=species_taxids,
        columns=two_odd_ids_sorted,
        dtype=int,
    )
    for tid in species_taxids:
        for two_odd_id in presence_by_tid.get(tid, set()):
            if two_odd_id in base.columns:
                base.at[tid, two_odd_id] = 1

    ncbi = NCBITaxa()

    species_counts_by_row: dict[int, int] | None = None

    if rank == "species":
        row_taxids = species_taxids
        # Species level is binary (0/1).
        matrix = base.astype(int)
        row_labels = ncbi.get_taxid_translator(row_taxids)
        # Don't add parentheses for species labels.
        species_counts_by_row = None
        ordered_rows = _topology_order(ncbi, row_taxids)
        matrix = matrix.loc[ordered_rows]
        rows_for_tree = ordered_rows
    else:
        # map species -> rank taxid
        sp_to_rank: dict[int, int] = {}
        for sp_tid in species_taxids:
            rt = _rank_taxid(ncbi, sp_tid, rank)
            if rt is not None:
                sp_to_rank[sp_tid] = rt

        if not sp_to_rank:
            log_line(log_path, f"[ANALYZE] Could not map any species to rank '{rank}'; skipping")
            return

        rank_taxids = sorted(set(sp_to_rank.values()))
        ordered_rows = _topology_order(ncbi, rank_taxids)

        # Count how many species are represented in each row taxon.
        species_counts_by_row = {}
        for sp_tid, rt in sp_to_rank.items():
            species_counts_by_row[rt] = int(species_counts_by_row.get(rt, 0)) + 1

        # aggregate to % of species within each rank having the 2ODD
        rows = []
        for rt in ordered_rows:
            sp_in_rank = [sp for sp, rt2 in sp_to_rank.items() if rt2 == rt]
            if not sp_in_rank:
                continue
            frac = base.loc[sp_in_rank].mean(axis=0)  # fraction of species with presence
            rows.append((rt, (frac * 100.0).astype(float)))

        matrix = pd.DataFrame({tid: row for tid, row in rows}).T
        matrix.index.name = rank
        matrix = matrix[two_odd_ids_sorted]
        row_labels = ncbi.get_taxid_translator(list(matrix.index))
        rows_for_tree = list(matrix.index)

    display_labels = _build_display_labels(row_labels, species_counts_by_row)

    # Build color mapping for the displayed tip labels.
    label_colors: dict[str, str] = {}
    present_groups: set[str] = set()
    for tid in rows_for_tree:
        group = _plant_group_for_taxid(ncbi, int(tid))
        if group:
            present_groups.add(group)
            color = GROUP_COLORS.get(group)
            if color:
                label = display_labels.get(int(tid))
                if label:
                    label_colors[label] = color

    # write outputs
    matrix_tsv = output_dir / f"presence_matrix_{rank}.tsv"
    matrix.to_csv(matrix_tsv, sep="\t", index=True)

    # Write a Newick tree for the row ordering
    tree_nwk = output_dir / f"presence_tree_{rank}.nwk"
    try:
        topo = ncbi.get_topology([int(x) for x in rows_for_tree], annotate=False)
        topo.write(outfile=str(tree_nwk))
        log_line(log_path, f"[ANALYZE] Wrote topology tree: {tree_nwk.name}")
    except Exception as e:
        log_line(log_path, f"[ANALYZE] Could not write topology tree: {e}")

    # plot (lazy import so pipeline can run without plotting deps unless enabled)
    try:
        import matplotlib.pyplot as plt
        from Bio import Phylo
    except Exception as e:
        log_line(
            log_path,
            f"[ANALYZE] Plotting skipped because matplotlib is unavailable: {e}. Matrix saved to {matrix_tsv}",
        )
        return

    # Build a tree with taxid leaf names so it matches the matrix order
    tree = None
    if tree_nwk.exists():
        try:
            tree = Phylo.read(str(tree_nwk), "newick")
            _make_ultrametric_biophylo(tree)
            _sort_biophylo_tree_by_terminal_order(tree, [int(x) for x in matrix.index])
        except Exception as e:
            log_line(log_path, f"[ANALYZE] Could not read {tree_nwk.name} with Bio.Phylo: {e}")
            tree = None

    n_rows = len(matrix)
    fig_height = max(10, 0.25 * n_rows)
    fig = plt.figure(figsize=(22, fig_height))
    gs = fig.add_gridspec(1, 2, width_ratios=[1, 3])

    ax_tree = fig.add_subplot(gs[0])
    if tree is not None:

        def _label(clade):
            if clade.name and str(clade.name).isdigit():
                return display_labels.get(int(clade.name), str(clade.name))
            return None

        Phylo.draw(
            tree,
            axes=ax_tree,
            do_show=False,
            label_func=_label,
            label_colors=label_colors,
        )
        ax_tree.set_axis_off()
    else:
        ax_tree.text(0.5, 0.5, "Tree unavailable", ha="center", va="center")
        ax_tree.set_axis_off()

    # Add a legend for plant-group colors in the bottom-left whitespace.
    try:
        from matplotlib.patches import Patch

        groups_in_plot = [g for g in GROUP_COLORS.keys() if g in present_groups]
        if groups_in_plot:
            handles = [Patch(facecolor=GROUP_COLORS[g], edgecolor="none", label=g) for g in groups_in_plot]
            fig.legend(
                handles=handles,
                title="Plant groups",
                loc="lower left",
                bbox_to_anchor=(0.02, 0.02),
                frameon=False,
                fontsize=8,
                title_fontsize=9,
                ncol=1,
                handlelength=1.0,
                handleheight=1.0,
            )
    except Exception:
        # Legend is cosmetic; never fail analysis if it can't be drawn.
        pass

    ax_heat = fig.add_subplot(gs[1])
    arr = matrix.to_numpy(dtype=float)

    # Heatmap styling differs by rank.
    is_species_rank = rank == "species"
    vmax = 1.0 if is_species_rank else 100.0
    im = ax_heat.imshow(
        arr,
        aspect="auto",
        interpolation="nearest",
        cmap="viridis",
        vmin=0.0,
        vmax=vmax,
    )
    ax_heat.set_xticks(range(len(matrix.columns)))
    ax_heat.set_xticklabels(
        _two_odd_xtick_labels(list(matrix.columns), major_char_info),
        rotation=90,
    )

    # Tree already shows row labels; avoid duplicate labels on heatmap.
    ax_heat.set_yticks([])
    ax_heat.tick_params(axis="y", left=False)

    cbar = fig.colorbar(im, ax=ax_heat, fraction=0.035, pad=0.02)
    cbar.set_label("Presence (0/1)" if is_species_rank else "% of species with 2ODD")

    if not is_species_rank:
        # Annotate values directly in cells when the matrix is reasonably sized.
        # This becomes extremely slow and unreadable for very large matrices.
        max_cells_to_annotate = int(analyze_cfg.get("heatmap_max_annotated_cells", 20000))
        n_cells = int(arr.shape[0] * arr.shape[1])
        if n_cells <= max_cells_to_annotate:
            for i in range(arr.shape[0]):
                for j in range(arr.shape[1]):
                    val = float(arr[i, j])
                    text = f"{val:.0f}"
                    # Choose contrasting text color based on background
                    txt_color = "white" if val >= 50.0 else "black"
                    ax_heat.text(j, i, text, ha="center", va="center", fontsize=6, color=txt_color)
        else:
            log_line(
                log_path,
                f"[ANALYZE] Skipping cell annotations: matrix has {n_cells} cells (limit {max_cells_to_annotate}).",
            )

    if is_species_rank:
        title = "2ODD presence by taxonomic species level\nx: (assoc. functions) 2ODD ID"
    else:
        title = (
            f"2ODD percentage presence by taxonomic {rank} level\n"
            "row (N)=#species; x: (assoc. functions) 2ODD ID"
        )
    ax_heat.set_title(title)

    plt.tight_layout()
    out_png = output_dir / f"presence_heatmap_{rank}.png"
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close(fig)

    log_line(log_path, f"[ANALYZE] Wrote {out_png.name} and {matrix_tsv.name}")
