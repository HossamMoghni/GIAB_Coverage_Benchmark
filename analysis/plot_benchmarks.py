#!/usr/bin/env python3
"""Regenerate every benchmarking figure from the committed hap.py summary tables.

The figures in results/figures/ are build artifacts, not hand-made images: this
script is the only thing that produces them, so a corrected benchmarking run
(see scripts/07_rerun_ontarget.sh) redraws the whole set with one command.

Usage
-----
    python analysis/plot_benchmarks.py
    python analysis/plot_benchmarks.py --results-dir results/happy_ontarget \\
                                       --out-dir results/figures_ontarget

Input is one hap.py ``*.summary.csv`` per depth, named ``<chrom>_<depth>x.summary.csv``.
Only the ``PASS`` filter rows are read, because HaplotypeCaller output here is
unfiltered and the ALL and PASS rows are identical; taking PASS keeps the script
correct if hard filtering is added later.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

# --- Design tokens ------------------------------------------------------------
# Categorical slots 1-3 of the validated default palette. Three is the cap for
# all-pairs colorblind separation, which is exactly what these charts need. The
# aqua slot sits below 3:1 against the surface, so every chart that uses it also
# carries visible direct labels (the relief rule).

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_MUTED = "#52514e"
GRID = "#e3e2de"

SERIES = ["#2a78d6", "#eb6834", "#1baf7a"]  # blue, orange, aqua
SEQUENTIAL = [  # single-hue blue ramp, light -> dark
    "#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec",
    "#5598e7", "#3987e5", "#2a78d6", "#256abf", "#1c5cab",
    "#184f95", "#104281", "#0d366b",
]

VARIANT_COLOR = {"SNP": SERIES[0], "INDEL": SERIES[1]}
METRIC_COLOR = {"Recall": SERIES[0], "Precision": SERIES[1], "F1-score": SERIES[2]}


def style() -> None:
    """Recessive axes, thin marks, text in ink tokens rather than series color."""
    plt.rcParams.update({
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        "axes.edgecolor": GRID,
        "axes.labelcolor": INK_MUTED,
        "axes.titlecolor": INK,
        "axes.titlesize": 13,
        "axes.titleweight": "bold",
        "axes.labelsize": 10,
        "axes.grid": True,
        "axes.axisbelow": True,
        "grid.color": GRID,
        "grid.linewidth": 0.8,
        "xtick.color": INK_MUTED,
        "ytick.color": INK_MUTED,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.frameon": False,
        "legend.fontsize": 9,
        "legend.labelcolor": INK_MUTED,
        "lines.linewidth": 2.0,
        "lines.markersize": 8,
        "font.size": 10,
    })


def despine(ax) -> None:
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.grid(axis="x", visible=False)


# --- Data loading -------------------------------------------------------------

def load(results_dir: Path) -> pd.DataFrame:
    """Read every per-depth summary CSV into one tidy frame.

    Returns columns: depth, variant_type, recall, precision, f1, tp, fp, fn, total.
    """
    files = sorted(results_dir.glob("summary/*.summary.csv"))
    if not files:
        sys.exit(f"no summary CSVs found under {results_dir}/summary/")

    rows = []
    for path in files:
        match = re.search(r"_(\d+)x\.summary\.csv$", path.name)
        if not match:
            print(f"  skipping unrecognised filename: {path.name}", file=sys.stderr)
            continue
        depth = int(match.group(1))

        table = pd.read_csv(path)
        table = table[table["Filter"] == "PASS"]
        for _, row in table.iterrows():
            rows.append({
                "depth": depth,
                "variant_type": row["Type"],
                "recall": row["METRIC.Recall"] * 100,
                "precision": row["METRIC.Precision"] * 100,
                "f1": row["METRIC.F1_Score"] * 100,
                "tp": int(row["TRUTH.TP"]),
                "fp": int(row["QUERY.FP"]),
                "fn": int(row["TRUTH.FN"]),
                "total": int(row["TRUTH.TOTAL"]),
            })

    frame = pd.DataFrame(rows).sort_values(["variant_type", "depth"]).reset_index(drop=True)
    print(f"  loaded {len(files)} depths: {sorted(frame['depth'].unique())}")
    return frame


# --- Individual figures -------------------------------------------------------

def fig_snp_performance(df: pd.DataFrame, out: Path) -> None:
    """Recall, precision and F1 for SNPs against coverage depth.

    Coverage is on a linear axis rather than evenly spaced categories, because
    the diminishing-returns claim in the report is a statement about real depth
    intervals and category spacing would fabricate it.
    """
    snp = df[df.variant_type == "SNP"]
    fig, ax = plt.subplots(figsize=(8, 5))

    for label, column in (("Recall", "recall"), ("Precision", "precision"), ("F1-score", "f1")):
        ax.plot(snp.depth, snp[column], marker="o", color=METRIC_COLOR[label],
                label=label, markeredgecolor=SURFACE, markeredgewidth=2, zorder=3)
        # Direct labels on every series end: identity is never color-alone, and
        # this is the relief the aqua slot's contrast WARN requires.
        last = snp.iloc[-1]
        ax.annotate(f"{label}  {last[column]:.1f}%",
                    xy=(last.depth, last[column]), xytext=(8, 0),
                    textcoords="offset points", va="center",
                    color=INK_MUTED, fontsize=9)

    ax.set_xlabel("Mean coverage depth (x)")
    ax.set_ylabel("Metric (%)")
    ax.set_title("SNP detection against coverage depth")
    ax.set_xlim(0, snp.depth.max() * 1.42)
    ax.set_ylim(0, 100)
    ax.set_xticks(sorted(snp.depth))
    ax.set_xticklabels([f"{d}x" for d in sorted(snp.depth)])
    ax.legend(loc="upper left")
    despine(ax)
    save(fig, out)


def fig_recall_snp_vs_indel(df: pd.DataFrame, out: Path) -> None:
    """Recall for SNPs and indels side by side, as grouped bars.

    Bars, not lines: the comparison people make here is per-depth SNP versus
    indel, and adjacent bars put the two values against a shared baseline.
    """
    depths = sorted(df.depth.unique())
    x = np.arange(len(depths))
    width = 0.38

    fig, ax = plt.subplots(figsize=(8, 5))
    for offset, vtype in ((-width / 2 - 0.01, "SNP"), (width / 2 + 0.01, "INDEL")):
        vals = [df[(df.depth == d) & (df.variant_type == vtype)].recall.iloc[0] for d in depths]
        bars = ax.bar(x + offset, vals, width, color=VARIANT_COLOR[vtype], label=vtype,
                      edgecolor=SURFACE, linewidth=2, zorder=3)
        ax.bar_label(bars, fmt="%.1f%%", padding=3, color=INK_MUTED, fontsize=9)

    ax.set_xticks(x, [f"{d}x" for d in depths])
    ax.set_xlabel("Mean coverage depth")
    ax.set_ylabel("Recall (%)")
    ax.set_title("Recall by variant class")
    ax.set_ylim(0, max(df.recall) * 1.25)
    ax.legend(loc="upper left")
    despine(ax)
    save(fig, out)


def fig_f1_snp_vs_indel(df: pd.DataFrame, out: Path) -> None:
    """F1 for SNPs and indels against depth."""
    fig, ax = plt.subplots(figsize=(8, 5))

    for vtype in ("SNP", "INDEL"):
        sub = df[df.variant_type == vtype]
        ax.plot(sub.depth, sub.f1, marker="o", color=VARIANT_COLOR[vtype], label=vtype,
                markeredgecolor=SURFACE, markeredgewidth=2, zorder=3)
        last = sub.iloc[-1]
        ax.annotate(f"{vtype}  {last.f1:.1f}%", xy=(last.depth, last.f1), xytext=(8, 0),
                    textcoords="offset points", va="center", color=INK_MUTED, fontsize=9)

    ax.set_xlabel("Mean coverage depth (x)")
    ax.set_ylabel("F1-score (%)")
    ax.set_title("F1-score by variant class")
    ax.set_xlim(0, df.depth.max() * 1.3)
    ax.set_ylim(0, max(df.f1) * 1.3)
    ax.set_xticks(sorted(df.depth), [f"{d}x" for d in sorted(df.depth)])
    ax.legend(loc="upper left")
    despine(ax)
    save(fig, out)


def fig_precision_recall(df: pd.DataFrame, out: Path) -> None:
    """Precision against recall, one trajectory per variant class.

    Each point is a depth, so the path shows the actual trade-off: recall moves a
    long way right while precision barely moves up.
    """
    fig, ax = plt.subplots(figsize=(8, 5))

    for vtype in ("SNP", "INDEL"):
        sub = df[df.variant_type == vtype]
        ax.plot(sub.recall, sub.precision, marker="o", color=VARIANT_COLOR[vtype],
                label=vtype, markeredgecolor=SURFACE, markeredgewidth=2, zorder=3)
        # The two trajectories converge near the origin, so the depth labels are
        # pushed to opposite sides of their markers to keep them legible.
        dy = 11 if vtype == "SNP" else -17
        for _, point in sub.iterrows():
            ax.annotate(f"{point.depth:.0f}x", xy=(point.recall, point.precision),
                        xytext=(0, dy), textcoords="offset points", ha="center",
                        color=INK_MUTED, fontsize=8)

    ax.set_xlabel("Recall (%)")
    ax.set_ylabel("Precision (%)")
    ax.set_title("Precision–recall trade-off across coverage depths")
    ax.set_xlim(0, max(df.recall) * 1.2)
    ax.set_ylim(0, 100)
    ax.legend(loc="lower right")
    despine(ax)
    save(fig, out)


def fig_coverage_efficiency(df: pd.DataFrame, out: Path) -> None:
    """True positives recovered per unit of coverage, as small multiples.

    SNP and indel counts differ by an order of magnitude. Two panels sharing the
    x axis keep both readable; a second y axis would make the two series look
    comparable when they are not.
    """
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))

    for ax, vtype in zip(axes, ("SNP", "INDEL")):
        sub = df[df.variant_type == vtype].sort_values("depth")
        efficiency = sub.tp.to_numpy() / sub.depth.to_numpy()
        bars = ax.bar([f"{d}x" for d in sub.depth], efficiency,
                      color=VARIANT_COLOR[vtype], width=0.6,
                      edgecolor=SURFACE, linewidth=2, zorder=3)
        ax.bar_label(bars, fmt="%.0f", padding=3, color=INK_MUTED, fontsize=9)
        ax.set_title(f"{vtype}s recovered per 1x of coverage")
        ax.set_xlabel("Mean coverage depth")
        ax.set_ylim(0, efficiency.max() * 1.25)
        despine(ax)

    axes[0].set_ylabel("True positives per 1x")
    fig.suptitle("Detection efficiency: marginal yield of each unit of depth",
                 color=INK, fontsize=13, fontweight="bold", y=1.0)
    # Two panels, one series each; the legend keeps the color-to-class mapping
    # explicit rather than leaving it to the panel titles alone.
    fig.legend(handles=[Line2D([], [], color=VARIANT_COLOR[v], lw=6, label=v)
                        for v in ("SNP", "INDEL")],
               loc="lower center", ncol=2, bbox_to_anchor=(0.5, -0.04))
    save(fig, out)


def fig_performance_heatmap(df: pd.DataFrame, out: Path) -> None:
    """All metrics at all depths in one grid.

    A single-hue sequential ramp encodes magnitude; every cell also carries its
    value, so the color is a scanning aid rather than the only channel.
    """
    depths = sorted(df.depth.unique())
    rows = [(v, m, c) for v in ("SNP", "INDEL")
            for m, c in (("Recall", "recall"), ("Precision", "precision"), ("F1-score", "f1"))]

    matrix = np.array([[df[(df.depth == d) & (df.variant_type == v)][col].iloc[0]
                        for d in depths] for v, _, col in rows])

    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("seq_blue", SEQUENTIAL)

    fig, ax = plt.subplots(figsize=(7.5, 5))
    mesh = ax.imshow(matrix, cmap=cmap, vmin=0, vmax=100, aspect="auto")

    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            value = matrix[i, j]
            # Flip the label to the surface color once the cell is dark enough
            # to swallow dark ink.
            ax.text(j, i, f"{value:.1f}", ha="center", va="center", fontsize=10,
                    color=SURFACE if value > 55 else INK)

    ax.set_xticks(range(len(depths)), [f"{d}x" for d in depths])
    ax.set_yticks(range(len(rows)), [f"{v}  {m}" for v, m, _ in rows])
    ax.set_title("Benchmarking metrics (%) across coverage depths")
    ax.grid(False)
    for side in ax.spines.values():
        side.set_visible(False)

    bar = fig.colorbar(mesh, ax=ax, pad=0.02)
    bar.set_label("Metric (%)", color=INK_MUTED)
    bar.outline.set_visible(False)
    save(fig, out)


def fig_dashboard(df: pd.DataFrame, out: Path) -> None:
    """One-page composite of the four charts the report leans on."""
    depths = sorted(df.depth.unique())
    fig, axes = plt.subplots(2, 2, figsize=(13, 9))

    # (A) SNP metrics against depth
    ax = axes[0][0]
    snp = df[df.variant_type == "SNP"]
    for label, column in (("Recall", "recall"), ("Precision", "precision"), ("F1-score", "f1")):
        ax.plot(snp.depth, snp[column], marker="o", color=METRIC_COLOR[label], label=label,
                markeredgecolor=SURFACE, markeredgewidth=2, zorder=3)
    ax.set_title("A · SNP metrics against depth")
    ax.set_xlabel("Coverage (x)")
    ax.set_ylabel("Metric (%)")
    ax.set_ylim(0, 100)
    ax.set_xticks(depths, [f"{d}x" for d in depths])
    ax.legend(loc="upper left")
    despine(ax)

    # (B) Recall by variant class
    ax = axes[0][1]
    x = np.arange(len(depths))
    width = 0.38
    for offset, vtype in ((-width / 2 - 0.01, "SNP"), (width / 2 + 0.01, "INDEL")):
        vals = [df[(df.depth == d) & (df.variant_type == vtype)].recall.iloc[0] for d in depths]
        bars = ax.bar(x + offset, vals, width, color=VARIANT_COLOR[vtype], label=vtype,
                      edgecolor=SURFACE, linewidth=2, zorder=3)
        ax.bar_label(bars, fmt="%.1f", padding=3, color=INK_MUTED, fontsize=8)
    ax.set_title("B · Recall by variant class")
    ax.set_xlabel("Coverage")
    ax.set_ylabel("Recall (%)")
    ax.set_xticks(x, [f"{d}x" for d in depths])
    ax.set_ylim(0, max(df.recall) * 1.3)
    ax.legend(loc="upper left")
    despine(ax)

    # (C) Precision-recall trajectory
    ax = axes[1][0]
    for vtype in ("SNP", "INDEL"):
        sub = df[df.variant_type == vtype]
        ax.plot(sub.recall, sub.precision, marker="o", color=VARIANT_COLOR[vtype], label=vtype,
                markeredgecolor=SURFACE, markeredgewidth=2, zorder=3)
        for _, point in sub.iterrows():
            ax.annotate(f"{point.depth:.0f}x", xy=(point.recall, point.precision),
                        xytext=(0, 10 if vtype == "SNP" else -17),
                        textcoords="offset points", ha="center",
                        color=INK_MUTED, fontsize=8)
    ax.set_title("C · Precision–recall trade-off")
    ax.set_xlabel("Recall (%)")
    ax.set_ylabel("Precision (%)")
    ax.set_xlim(0, max(df.recall) * 1.2)
    ax.set_ylim(0, 100)
    ax.legend(loc="lower right")
    despine(ax)

    # (D) Variant counts
    ax = axes[1][1]
    table = df[df.variant_type == "SNP"].sort_values("depth")
    x = np.arange(len(depths))
    width = 0.27
    for i, (label, column) in enumerate((("True positives", "tp"), ("False positives", "fp"))):
        bars = ax.bar(x + (i - 0.5) * (width + 0.02), table[column], width,
                      color=SERIES[i], label=label, edgecolor=SURFACE, linewidth=2, zorder=3)
        ax.bar_label(bars, fmt="%d", padding=3, color=INK_MUTED, fontsize=8)
    ax.set_title("D · SNP call counts")
    ax.set_xlabel("Coverage")
    ax.set_ylabel("Variants")
    ax.set_xticks(x, [f"{d}x" for d in depths])
    ax.set_ylim(0, table.tp.max() * 1.25)
    ax.legend(loc="upper left")
    despine(ax)

    fig.suptitle("Variant calling performance across coverage depths · HG002 exome, chr22",
                 color=INK, fontsize=15, fontweight="bold")
    save(fig, out)


# --- Plumbing -----------------------------------------------------------------

def save(fig, path: Path) -> None:
    fig.tight_layout()
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {path}")


FIGURES = {
    "snp_performance_vs_coverage.png": fig_snp_performance,
    "recall_snp_vs_indel.png": fig_recall_snp_vs_indel,
    "f1_snp_vs_indel.png": fig_f1_snp_vs_indel,
    "precision_recall.png": fig_precision_recall,
    "coverage_efficiency.png": fig_coverage_efficiency,
    "performance_heatmap.png": fig_performance_heatmap,
    "dashboard.png": fig_dashboard,
}


def main() -> None:
    repo = Path(__file__).resolve().parent.parent
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--results-dir", type=Path, default=repo / "results" / "happy",
                        help="directory holding summary/*.summary.csv")
    parser.add_argument("--out-dir", type=Path, default=repo / "results" / "figures",
                        help="where to write the PNGs")
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    style()

    print(f"Reading {args.results_dir}")
    df = load(args.results_dir)

    print(f"Writing figures to {args.out_dir}")
    for name, builder in FIGURES.items():
        builder(df, args.out_dir / name)

    # A machine-readable roll-up, handy for the README table and for anyone who
    # wants the numbers without parsing hap.py's 60-column output.
    tidy = args.out_dir.parent / "benchmark_summary.csv"
    df.to_csv(tidy, index=False)
    print(f"  wrote {tidy}")


if __name__ == "__main__":
    main()
