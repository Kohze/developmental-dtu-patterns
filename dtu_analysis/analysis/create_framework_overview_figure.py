#!/usr/bin/env python3
"""Create the general-framework and mouse-brain application overview figure."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch


ROOT = Path(__file__).resolve().parents[2]
DTU_ROOT = ROOT / "dtu_analysis"
OUTPUTS = (
    ROOT / "figures" / "figure01_framework_overview",
    DTU_ROOT / "figures" / "figure01_framework_overview",
)

COLORS = {
    "ink": "#263746",
    "muted": "#657786",
    "grid": "#DCE3E8",
    "orange": "#D55E00",
    "blue": "#0072B2",
    "green": "#009E73",
    "light_orange": "#FBEDE4",
    "light_blue": "#E8F3F8",
    "light_green": "#E5F5EF",
    "light_gray": "#F2F5F7",
}


def rounded_box(ax, xy, width, height, text, facecolor, edgecolor, fontsize=10.5):
    box = FancyBboxPatch(
        xy,
        width,
        height,
        boxstyle="round,pad=0.018,rounding_size=0.025",
        linewidth=1.8,
        edgecolor=edgecolor,
        facecolor=facecolor,
        transform=ax.transAxes,
        clip_on=False,
    )
    ax.add_patch(box)
    ax.text(
        xy[0] + width / 2,
        xy[1] + height / 2,
        text,
        transform=ax.transAxes,
        ha="center",
        va="center",
        color=COLORS["ink"],
        fontsize=fontsize,
        linespacing=1.3,
    )


def add_arrow(ax, start, end):
    ax.annotate(
        "",
        xy=end,
        xytext=start,
        xycoords=ax.transAxes,
        textcoords=ax.transAxes,
        arrowprops=dict(arrowstyle="-|>", color=COLORS["muted"], lw=1.8),
    )


def load_results():
    sensitivity = pd.read_csv(DTU_ROOT / "tables" / "transient_regional_sensitivity.csv")
    dependence = pd.read_csv(
        DTU_ROOT / "tables" / "transient_regional_dependence_sensitivity.csv"
    )
    ranking = pd.read_csv(DTU_ROOT / "tables" / "transient_regional_top_candidates.csv")
    primary = dependence.loc[dependence["method"] == "primary"].iloc[0]
    robust = dependence.loc[dependence["method"] != "primary"].iloc[0]
    return sensitivity, primary, robust, ranking


def create_figure():
    sensitivity, primary, robust, ranking = load_results()

    fig = plt.figure(figsize=(15.6, 10.2), facecolor="white")
    grid = fig.add_gridspec(
        2,
        2,
        height_ratios=[1.02, 1.0],
        width_ratios=[1.1, 0.9],
        hspace=0.32,
        wspace=0.24,
    )

    # Panel A: the general pattern, defined independently of this dataset.
    ax_a = fig.add_subplot(grid[0, :])
    stages = np.arange(5)
    target = np.array([0.48, 0.50, 0.84, 0.49, 0.50])
    reference_1 = np.array([0.47, 0.49, 0.50, 0.48, 0.51])
    reference_2 = np.array([0.49, 0.48, 0.48, 0.50, 0.49])
    ax_a.axvspan(1.72, 2.28, color="#FFF6C7", zorder=0)
    ax_a.plot(stages, target, "-o", color=COLORS["orange"], lw=3.0, ms=7,
              label="Focal group")
    ax_a.plot(stages, reference_1, "--o", color=COLORS["blue"], lw=2.4, ms=6,
              label="Comparison group 1")
    ax_a.plot(stages, reference_2, ":o", color=COLORS["green"], lw=2.6, ms=6,
              label="Comparison group 2")
    ax_a.set_xlim(-0.15, 4.15)
    ax_a.set_ylim(0.35, 0.92)
    ax_a.set_xticks(stages, ["flank -2", "flank -1", "candidate stage", "flank +1", "flank +2"])
    ax_a.set_ylabel("Relative transcript fraction", fontsize=11)
    ax_a.grid(axis="y", color=COLORS["grid"], lw=1)
    ax_a.spines[["top", "right"]].set_visible(False)
    ax_a.tick_params(colors=COLORS["muted"])
    ax_a.set_title(
        "A  General target: a group-specific divergence bounded by temporal reconvergence",
        loc="left",
        fontsize=16,
        fontweight="bold",
        color=COLORS["ink"],
        pad=13,
    )
    ax_a.text(
        2,
        0.885,
        "same-direction separation\nfrom both comparison groups",
        ha="center",
        va="top",
        fontsize=10.5,
        color=COLORS["ink"],
    )
    ax_a.text(
        0.5,
        0.385,
        "comparison groups agree",
        ha="center",
        fontsize=9.5,
        color=COLORS["muted"],
    )
    ax_a.text(
        3.5,
        0.385,
        "immediate flanks reconverge",
        ha="center",
        fontsize=9.5,
        color=COLORS["muted"],
    )
    ax_a.legend(
        frameon=False,
        ncol=3,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.26),
        fontsize=10,
    )

    # Panel B: model-agnostic post-inference decision layer.
    ax_b = fig.add_subplot(grid[1, 0])
    ax_b.axis("off")
    ax_b.set_title(
        "B  Model-agnostic post-inference decision layer",
        loc="left",
        fontsize=16,
        fontweight="bold",
        color=COLORS["ink"],
        pad=10,
    )
    y = 0.57
    width = 0.168
    height = 0.25
    xs = [0.00, 0.205, 0.410, 0.615, 0.820]
    labels = [
        "Replicate-level\ntranscript fractions\nand ordered stages",
        "Stage-specific\ncomponent tests\nfrom a DTU engine",
        "Direction, effect,\ncomparison-group\nagreement",
        "Immediate-flank\nreconvergence and\nreplicate separation",
        "Candidate episodes\nand programmatically\nranked panel",
    ]
    faces = [
        COLORS["light_gray"],
        COLORS["light_blue"],
        COLORS["light_orange"],
        COLORS["light_green"],
        COLORS["light_orange"],
    ]
    edges = [COLORS["muted"], COLORS["blue"], COLORS["orange"], COLORS["green"], COLORS["orange"]]
    for x, label, face, edge in zip(xs, labels, faces, edges):
        rounded_box(ax_b, (x, y), width, height, label, face, edge, fontsize=9.4)
    for left, right in zip(xs[:-1], xs[1:]):
        add_arrow(ax_b, (left + width + 0.006, y + height / 2), (right - 0.008, y + height / 2))
    ax_b.text(
        0.0,
        0.35,
        "Required design: at least three ordered stages, one focal group,\n"
        "at least two comparison groups and biological replication.",
        transform=ax_b.transAxes,
        fontsize=10.5,
        color=COLORS["ink"],
        va="top",
    )
    ax_b.text(
        0.0,
        0.18,
        "The framework organises statistically supported components into\n"
        "candidate temporal episodes; it does not replace the underlying DTU model.",
        transform=ax_b.transAxes,
        fontsize=10.5,
        color=COLORS["muted"],
        va="top",
    )

    # Panel C: application-specific outcomes and sensitivity.
    ax_c = fig.add_subplot(grid[1, 1])
    ax_c.axis("off")
    ax_c.set_title(
        "C  Application to the mouse developmental brain archive",
        loc="left",
        fontsize=16,
        fontweight="bold",
        color=COLORS["ink"],
        pad=10,
    )
    rounded_box(
        ax_c,
        (0.01, 0.68),
        0.45,
        0.20,
        "12,517 isoforms\n4,577 multi-isoform genes",
        COLORS["light_blue"],
        COLORS["blue"],
        fontsize=11.5,
    )
    rounded_box(
        ax_c,
        (0.53, 0.68),
        0.45,
        0.20,
        f"{int(primary['replicate_separated_episodes']):,} candidate episodes\n"
        f"{int(primary['genes']):,} genes",
        COLORS["light_orange"],
        COLORS["orange"],
        fontsize=11.5,
    )
    add_arrow(ax_c, (0.465, 0.78), (0.525, 0.78))
    rounded_box(
        ax_c,
        (0.01, 0.40),
        0.45,
        0.20,
        f"{int(robust['replicate_separated_episodes']):,} episodes / "
        f"{int(robust['genes']):,} genes\nunder BY/Bonferroni sensitivity",
        COLORS["light_green"],
        COLORS["green"],
        fontsize=10.8,
    )
    rounded_box(
        ax_c,
        (0.53, 0.40),
        0.45,
        0.20,
        "Six programmatically\nhighest-ranked reciprocal genes",
        COLORS["light_orange"],
        COLORS["orange"],
        fontsize=10.8,
    )
    ax_c.text(
        0.01,
        0.29,
        "Leading genes: " + ", ".join(ranking["gene_name"].tolist()),
        transform=ax_c.transAxes,
        fontsize=10.5,
        color=COLORS["ink"],
    )
    sensitivity_range = (
        int(sensitivity["episodes"].min()),
        int(sensitivity["episodes"].max()),
    )
    ax_c.text(
        0.01,
        0.17,
        f"One-at-a-time threshold sensitivities: {sensitivity_range[0]:,}--"
        f"{sensitivity_range[1]:,} episodes.",
        transform=ax_c.transAxes,
        fontsize=10.5,
        color=COLORS["muted"],
    )
    ax_c.text(
        0.01,
        0.06,
        "All six leading genes were unchanged under the dependence-robust analysis.",
        transform=ax_c.transAxes,
        fontsize=10.5,
        color=COLORS["muted"],
    )

    fig.subplots_adjust(left=0.065, right=0.975, top=0.95, bottom=0.075)
    for stem in OUTPUTS:
        stem.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
        fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    create_figure()
