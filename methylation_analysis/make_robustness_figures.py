"""Publication figures for the replicate-aware methylation--isoform audit."""

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle
import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
FIGURES = HERE / "figures"
FIGURES.mkdir(exist_ok=True)

STAGES = ["10.5", "11.5", "12.5", "13.5", "14.5", "15.5", "16.5", "0"]
X = np.arange(len(STAGES))
STAGE_LABELS = [f"E{x}" if x != "0" else "P0" for x in STAGES]

NAVY = "#16324F"
BLUE = "#2878B5"
CYAN = "#50A6C2"
ORANGE = "#E07A3F"
RED = "#B5433C"
GOLD = "#D5A021"
GREY = "#89939E"
LIGHT = "#EEF2F5"
DARK = "#263238"

mpl.rcParams.update(
    {
        "font.family": "Arial",
        "font.size": 12.5,
        "axes.titlesize": 13.5,
        "axes.labelsize": 12.5,
        "xtick.labelsize": 11,
        "ytick.labelsize": 11,
        "legend.fontsize": 10.5,
        "axes.titleweight": "bold",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "savefig.bbox": "tight",
    }
)


def panel(ax, letter):
    ax.text(
        -0.12,
        1.08,
        letter,
        transform=ax.transAxes,
        fontsize=16,
        fontweight="bold",
        va="top",
    )


def save(fig, name):
    fig.savefig(FIGURES / f"{name}.pdf", dpi=300)
    fig.savefig(FIGURES / f"{name}.png", dpi=300)
    plt.close(fig)


def stage_profile(ax, frame, value, color, ylabel, title, ylim=None):
    for rep, marker in [(1, "o"), (2, "s")]:
        part = frame[frame["replicate"] == rep].set_index("stage").reindex(STAGES)
        ax.scatter(
            X,
            part[value],
            s=28,
            marker=marker,
            facecolor="white",
            edgecolor=color,
            linewidth=1.2,
            zorder=3,
            label=f"replicate {rep}",
        )
    mean = frame.groupby("stage")[value].mean().reindex(STAGES)
    ax.plot(X, mean, color=color, linewidth=2.2, zorder=2)
    ax.set_xticks(X, STAGE_LABELS, rotation=35, ha="right")
    ax.set_ylabel(ylabel)
    ax.set_title(title, loc="left")
    ax.grid(axis="y", color="#D9E0E5", linewidth=0.6)
    if ylim is not None:
        ax.set_ylim(*ylim)


def figure_audit(audit):
    fig = plt.figure(figsize=(11, 8.5))
    gs = fig.add_gridspec(2, 2, height_ratios=[0.78, 1.22], hspace=0.52, wspace=0.46)
    ax_flow = fig.add_subplot(gs[0, :])
    ax_scatter = fig.add_subplot(gs[1, 0])
    ax_sensitivity = fig.add_subplot(gs[1, 1])

    ax_flow.axis("off")
    panel(ax_flow, "a")
    boxes = [
        (0.02, "11,002", "archived tissue–region–\nisoform tests", NAVY),
        (0.365, "131", "region-wise BH\n(82 genes)", ORANGE),
        (0.71, "3", "parametric global BH\n(2 genes; no overlap)", BLUE),
    ]
    for xpos, number, text, color in boxes:
        patch = FancyBboxPatch(
            (xpos, 0.22),
            0.27,
            0.58,
            boxstyle="round,pad=0.016,rounding_size=0.02",
            facecolor="white",
            edgecolor=color,
            linewidth=2,
            transform=ax_flow.transAxes,
        )
        ax_flow.add_patch(patch)
        ax_flow.text(
            xpos + 0.135,
            0.58,
            number,
            ha="center",
            va="center",
            transform=ax_flow.transAxes,
            fontsize=24,
            fontweight="bold",
            color=color,
        )
        ax_flow.text(
            xpos + 0.135,
            0.34,
            text,
            ha="center",
            va="center",
            transform=ax_flow.transAxes,
            fontsize=11.5,
        )
    for x0, x1 in [(0.29, 0.365), (0.635, 0.71)]:
        ax_flow.add_patch(
            FancyArrowPatch(
                (x0, 0.51),
                (x1, 0.51),
                arrowstyle="-|>",
                mutation_scale=13,
                linewidth=1.5,
                color=GREY,
                transform=ax_flow.transAxes,
            )
        )
    ax_flow.text(
        0.5,
        0.04,
        "The parametric sets are disjoint; stage randomisation selects no row.",
        ha="center",
        transform=ax_flow.transAxes,
        color=DARK,
        fontweight="bold",
    )

    archived = audit["archived_q"] < 0.05
    reconstructed = audit["raw_global_q"] < 0.05
    ax_scatter.scatter(
        audit.loc[~archived & ~reconstructed, "archived_r"],
        audit.loc[~archived & ~reconstructed, "raw_weighted_r"],
        s=8,
        alpha=0.18,
        color=GREY,
        linewidth=0,
        rasterized=True,
    )
    ax_scatter.scatter(
        audit.loc[archived, "archived_r"],
        audit.loc[archived, "raw_weighted_r"],
        s=18,
        alpha=0.75,
        color=ORANGE,
        label="archived FDR",
    )
    ax_scatter.scatter(
        audit.loc[reconstructed, "archived_r"],
        audit.loc[reconstructed, "raw_weighted_r"],
        s=42,
        color=BLUE,
        edgecolor="white",
        linewidth=0.7,
        label="reconstructed FDR",
        zorder=4,
    )
    ax_scatter.axline((-1, -1), (1, 1), color=DARK, linewidth=0.8, linestyle="--")
    ax_scatter.axhline(0, color="#C9D0D5", linewidth=0.7)
    ax_scatter.axvline(0, color="#C9D0D5", linewidth=0.7)
    ax_scatter.set(xlim=(-1.04, 1.04), ylim=(-1.04, 1.04))
    ax_scatter.set_xlabel("Archived Pearson correlation")
    ax_scatter.set_ylabel("Reconstructed correlation")
    ax_scatter.set_title("Effect estimates are not\nreconstruction-invariant", loc="left")
    ax_scatter.legend(frameon=False, loc="lower right", fontsize=11)
    for gene, dx, dy in [("Gnao1", -0.02, 0.08), ("Taok3", 0.02, -0.12)]:
        point = audit[reconstructed & (audit["gene_id"] == gene)].iloc[0]
        ax_scatter.annotate(
            gene,
            (point.archived_r, point.raw_weighted_r),
            xytext=(point.archived_r + dx, point.raw_weighted_r + dy),
            fontsize=11,
            fontstyle="italic",
            arrowprops={"arrowstyle": "-", "color": GREY, "lw": 0.7},
        )
    panel(ax_scatter, "b")

    counts = {
        "Archived\n(region BH)": int(archived.sum()),
        "Parametric\n(global)": int(reconstructed.sum()),
        "Stage randomisation\n(global)": int(
            (audit["raw_weighted_exact_global_q"] < 0.05).sum()
        ),
        "Unweighted\n(global)": int((audit["unweighted_global_q"] < 0.05).sum()),
        "Strand-aware\n(global)": int(
            (audit["strand_sensitive_global_q"] < 0.05).sum()
        ),
    }
    positions = np.arange(len(counts))
    bars = ax_sensitivity.barh(
        positions,
        list(counts.values()),
        color=[ORANGE, BLUE, RED, CYAN, GREY],
        height=0.62,
    )
    ax_sensitivity.set_yticks(positions, counts.keys())
    ax_sensitivity.tick_params(axis="y", labelsize=10)
    ax_sensitivity.invert_yaxis()
    ax_sensitivity.set_xlabel("Rows at FDR < 0.05")
    ax_sensitivity.set_title("Discoveries depend on the\nanalysis specification", loc="left")
    ax_sensitivity.set_xlim(0, 145)
    ax_sensitivity.grid(axis="x", color="#D9E0E5", linewidth=0.6)
    for bar, value in zip(bars, counts.values()):
        ax_sensitivity.text(
            max(value, 0) + 3,
            bar.get_y() + bar.get_height() / 2,
            str(value),
            va="center",
            fontweight="bold",
        )
    ax_sensitivity.text(
        82,
        3.35,
        "No row survives randomisation FDR\nor all reasonable summaries",
        ha="center",
        va="center",
        color=RED,
        fontweight="bold",
    )
    panel(ax_sensitivity, "c")
    fig.suptitle(
        "Replicate-aware reconstruction replaces—not confirms—the archived discovery set",
        x=0.07,
        ha="left",
        fontsize=16.5,
        fontweight="bold",
        color=NAVY,
    )
    save(fig, "figure1_inferential_audit")


def figure_gnao1(methylation, usage):
    fig, axes = plt.subplots(2, 2, figsize=(11, 9.8))
    fig.subplots_adjust(hspace=0.70, wspace=0.48, top=0.88)
    m = methylation[
        (methylation.gene_id == "Gnao1")
        & (methylation.region == "upstream")
        & (methylation.tissue == "ForeBrain")
    ]
    u = usage[
        usage.isoform_id.isin(["NM_001113384", "NM_010308"])
        & (usage.tissue == "ForeBrain")
    ]
    stage_profile(
        axes[0, 0],
        m,
        "mean_weighted",
        BLUE,
        "Weighted CpG methylation (%)",
        "Forebrain 6-kb upstream methylation\nfalls during development",
        (0, 7.2),
    )
    axes[0, 0].legend(frameon=False, fontsize=10.5, ncol=2)
    panel(axes[0, 0], "a")

    for isoform, color, label in [
        ("NM_001113384", ORANGE, "NM_001113384"),
        ("NM_010308", BLUE, "NM_010308"),
    ]:
        part = u[u.isoform_id == isoform]
        mean = part.groupby("stage").rep_if.mean().reindex(STAGES)
        axes[0, 1].plot(X, mean, marker="o", color=color, lw=2.2, label=label)
        for rep in [1, 2]:
            points = part[part.replicate == rep].set_index("stage").reindex(STAGES)
            axes[0, 1].scatter(X, points.rep_if, s=18, facecolor="white", edgecolor=color)
    axes[0, 1].set_xticks(X, STAGE_LABELS, rotation=35, ha="right")
    axes[0, 1].set_ylabel("Relative isoform fraction")
    axes[0, 1].set_ylim(-0.03, 1.04)
    axes[0, 1].set_title(
        "Composition shifts from\nNM_001113384 to NM_010308", loc="left"
    )
    axes[0, 1].legend(frameon=False)
    axes[0, 1].grid(axis="y", color="#D9E0E5", linewidth=0.6)
    panel(axes[0, 1], "b")

    for isoform, color, label in [
        ("NM_001113384", ORANGE, "NM_001113384"),
        ("NM_010308", BLUE, "NM_010308"),
    ]:
        part = u[u.isoform_id == isoform]
        mean = part.groupby("stage").expression.mean().reindex(STAGES)
        axes[1, 0].plot(X, mean, marker="o", color=color, lw=2.2, label=label)
        for rep in [1, 2]:
            points = part[part.replicate == rep].set_index("stage").reindex(STAGES)
            axes[1, 0].scatter(
                X, points.expression, s=18, facecolor="white", edgecolor=color
            )
    axes[1, 0].set_xticks(X, STAGE_LABELS, rotation=35, ha="right")
    axes[1, 0].set_ylabel("Archived transcript expression")
    axes[1, 0].set_title("NM_010308 drives the rise\nin absolute abundance", loc="left")
    axes[1, 0].legend(frameon=False, loc="center right")
    axes[1, 0].grid(axis="y", color="#D9E0E5", linewidth=0.6)
    panel(axes[1, 0], "c")

    ax = axes[1, 1]
    ax.axis("off")
    ax.set_title("Evidence hierarchy for the\nGnao1 association", loc="left")
    y_positions = [0.78, 0.49, 0.20]
    items = [
        (
            "Observed",
            "6-kb upstream methylation ↓;\ntotal Gnao1 and NM_010308 ↑",
            BLUE,
        ),
        (
            "Reproduced pattern",
            "NM_001113384-to-NM_010308\nredistribution across datasets",
            GOLD,
        ),
        (
            "Mechanistic hypotheses",
            "local methylation, transcript control\nand cell composition",
            RED,
        ),
    ]
    for y, (heading, text, color) in zip(y_positions, items):
        ax.add_patch(
            FancyBboxPatch(
                (0.03, y - 0.105),
                0.92,
                0.21,
                boxstyle="round,pad=0.012",
                facecolor="white",
                edgecolor=color,
                linewidth=1.5,
                transform=ax.transAxes,
            )
        )
        ax.text(0.06, y + 0.050, heading, transform=ax.transAxes, color=color, fontweight="bold")
        ax.text(
            0.06,
            y - 0.035,
            text,
            transform=ax.transAxes,
            va="center",
            fontsize=10.2,
            linespacing=1.15,
        )
    for y0, y1 in zip(y_positions[:-1], y_positions[1:]):
        ax.add_patch(
            FancyArrowPatch(
                (0.49, y0 - 0.115),
                (0.49, y1 + 0.115),
                arrowstyle="-|>",
                mutation_scale=11,
                color=GREY,
                transform=ax.transAxes,
            )
        )
    panel(ax, "d")
    fig.suptitle(
        "Gnao1 links developmental transcript redistribution to methylation context",
        x=0.08,
        ha="left",
        fontsize=16.5,
        fontweight="bold",
        color=NAVY,
    )
    save(fig, "figure2_gnao1_composition")


def figure_taok3(methylation, usage, targeted):
    fig = plt.figure(figsize=(11, 9.2))
    gs = fig.add_gridspec(2, 2, hspace=0.62, wspace=0.48)
    axes = np.array(
        [
            [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])],
            [fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])],
        ]
    )
    m = methylation[
        (methylation.gene_id == "Taok3")
        & (methylation.region == "gene_body")
        & (methylation.tissue == "ForeBrain")
    ]
    u = usage[
        (usage.isoform_id == "NM_001199685") & (usage.tissue == "ForeBrain")
    ]
    stage_profile(
        axes[0, 0],
        m,
        "mean_weighted",
        BLUE,
        "Weighted CpG methylation (%)",
        "Taok3 gene-body methylation",
        (68, 79),
    )
    axes[0, 0].legend(frameon=False, fontsize=10.5, ncol=2)
    panel(axes[0, 0], "a")
    stage_profile(
        axes[0, 1],
        u,
        "rep_if",
        ORANGE,
        "Relative isoform fraction",
        "NM_001199685 relative usage",
        (0.3, 1.0),
    )
    panel(axes[0, 1], "b")

    row = targeted[
        (targeted.gene_id == "Taok3") & (targeted.isoform_id == "NM_001199685")
    ].iloc[0]
    labels = [
        "coverage-weighted",
        "unweighted CpGs",
        "linear detrending",
        "first differences",
        "replicate resampling",
    ]
    estimates = [
        row.raw_weighted_r,
        row.raw_unweighted_r,
        row.linear_detrended_r,
        row.first_difference_r,
        row.bootstrap_median_r,
    ]
    low = [np.nan, np.nan, np.nan, np.nan, row.bootstrap_low_r]
    high = [np.nan, np.nan, np.nan, np.nan, row.bootstrap_high_r]
    ypos = np.arange(len(labels))[::-1]
    axes[1, 0].axvline(0, color=DARK, linewidth=0.8)
    axes[1, 0].scatter(estimates, ypos, color=[BLUE, CYAN, GOLD, ORANGE, NAVY], s=45, zorder=3)
    axes[1, 0].hlines(
        ypos[-1],
        low[-1],
        high[-1],
        color=NAVY,
        linewidth=2,
    )
    axes[1, 0].set_yticks(ypos, labels)
    axes[1, 0].set_xlim(-0.05, 1.05)
    axes[1, 0].set_xlabel("Correlation estimate")
    axes[1, 0].set_title("Temporal association survives\nseveral stress tests", loc="left")
    axes[1, 0].grid(axis="x", color="#D9E0E5", linewidth=0.6)
    panel(axes[1, 0], "c")

    ax = axes[1, 1]
    ax.axis("off")
    ax.set_title("Why Taok3 remains a hypothesis,\nnot a mechanism", loc="left")
    items = [
        ("Strength", "detrended r = 0.968; first-difference r = 0.971", BLUE),
        ("Fragility", "weighted q = 0.038, but unweighted global q > 0.05", ORANGE),
        ("Replication", "WGBS profile reliability r = 0.343", RED),
        ("Annotation", "validated transcripts encode the same protein;\na methylation–start-site route is untested", GREY),
    ]
    for i, (head, text, color) in enumerate(items):
        y = 0.86 - i * 0.22
        ax.add_patch(Rectangle((0.04, y - 0.06), 0.025, 0.12, color=color, transform=ax.transAxes))
        ax.text(0.09, y + 0.026, head, color=color, fontweight="bold", transform=ax.transAxes)
        ax.text(0.09, y - 0.035, text, fontsize=10.5, transform=ax.transAxes, va="center")
    panel(ax, "d")
    fig.suptitle(
        "Taok3 is temporally coherent but sensitive to CpG weighting and replicate noise",
        x=0.08,
        ha="left",
        fontsize=16.5,
        fontweight="bold",
        color=NAVY,
    )
    save(fig, "figure3_taok3_sensitivity")


def figure_archived_robustness(targeted):
    archived = targeted[targeted["archived_q"] < 0.05].copy()
    fig, axes = plt.subplots(1, 3, figsize=(12, 5.6))
    fig.subplots_adjust(wspace=0.58, top=0.76)

    labels = [
        "Archived\nFDR",
        "Sign\nretained",
        "Nominal\np<0.05",
        "Stage shuffle\np<0.05",
    ]
    values = [
        len(archived),
        int(archived.archived_sign_reproduced.sum()),
        int(archived.raw_nominal_reproduced.sum()),
        int((archived.stage_permutation_p < 0.05).sum()),
    ]
    colors = [ORANGE, GOLD, CYAN, BLUE]
    bars = axes[0].bar(range(4), values, color=colors, width=0.7)
    axes[0].set_xticks(range(4), labels, fontsize=10)
    axes[0].set_ylabel("Archived rows")
    axes[0].set_ylim(0, 145)
    axes[0].set_title("Selection attrition\nunder reanalysis", loc="left")
    axes[0].grid(axis="y", color="#D9E0E5", linewidth=0.6)
    for bar, value in zip(bars, values):
        axes[0].text(bar.get_x() + bar.get_width() / 2, value + 3, value, ha="center", fontweight="bold")
    panel(axes[0], "a")

    axes[1].scatter(
        archived.raw_weighted_r,
        archived.linear_detrended_r,
        s=20,
        alpha=0.65,
        color=GOLD,
        label="linear detrending",
    )
    axes[1].scatter(
        archived.raw_weighted_r,
        archived.first_difference_r,
        s=20,
        alpha=0.65,
        color=BLUE,
        label="first differences",
    )
    axes[1].axline((-1, -1), (1, 1), color=DARK, linestyle="--", linewidth=0.8)
    axes[1].axhline(0, color="#C9D0D5", linewidth=0.7)
    axes[1].axvline(0, color="#C9D0D5", linewidth=0.7)
    axes[1].set(xlim=(-1.04, 1.04), ylim=(-1.04, 1.04))
    axes[1].set_xlabel("Reconstructed raw correlation")
    axes[1].set_ylabel("Temporal sensitivity correlation")
    axes[1].set_title("Stage trends explain\nmany large effects", loc="left")
    axes[1].legend(frameon=False, fontsize=10.5)
    panel(axes[1], "b")

    sc = axes[2].scatter(
        archived.methylation_replicate_reliability,
        archived.usage_replicate_reliability,
        c=archived.hostile_support_score,
        cmap="viridis",
        s=28,
        alpha=0.78,
        vmin=0,
        vmax=10,
    )
    axes[2].axhline(0.5, color=GREY, linestyle=":", linewidth=1)
    axes[2].axvline(0.5, color=GREY, linestyle=":", linewidth=1)
    axes[2].set(xlim=(-1.03, 1.03), ylim=(-1.03, 1.03))
    axes[2].set_xlabel("WGBS replicate-profile reliability")
    axes[2].set_ylabel("RNA replicate-profile reliability")
    axes[2].set_title("Reliability is uneven\nacross both assays", loc="left", pad=10)
    cbar = fig.colorbar(sc, ax=axes[2], fraction=0.046, pad=0.04)
    cbar.set_label("diagnostic support score")
    axes[2].text(
        -0.18,
        1.08,
        "c",
        transform=axes[2].transAxes,
        fontsize=16,
        fontweight="bold",
        va="top",
    )

    fig.suptitle(
        "The original 131 rows retain patterns, but not confirmatory multiplicity control",
        x=0.06,
        ha="left",
        fontsize=16.5,
        fontweight="bold",
        color=NAVY,
    )
    save(fig, "figure4_archived_robustness")


def main():
    audit = pd.read_csv(RESULTS / "replicate_robustness_all.csv")
    targeted = pd.read_csv(RESULTS / "replicate_robustness_archived_hits.csv")
    methylation = pd.read_csv(
        RESULTS / "replicate_level_methylation.csv", dtype={"stage": str}
    )
    usage = pd.read_csv(
        RESULTS / "replicate_level_isoform_usage.csv", dtype={"stage": str}
    )
    figure_audit(audit)
    figure_gnao1(methylation, usage)
    figure_taok3(methylation, usage, targeted)
    figure_archived_robustness(targeted)
    print("Wrote four audit figures to", FIGURES)


if __name__ == "__main__":
    main()
