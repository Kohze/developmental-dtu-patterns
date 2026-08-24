"""Generate manuscript figures from the archived thesis result tables.

This script is intentionally descriptive. It does not refit the underlying
WGBS or transcript-usage models and therefore does not convert the archived
correlations into causal or replicate-aware evidence.
"""

import argparse
import os
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT_DIR = Path(os.environ.get("METHYLATION_INPUT_DIR", ROOT)).expanduser()
OUT = Path(__file__).resolve().parent / "figures"
OUT.mkdir(exist_ok=True)

COLORS = {
    "blue": "#0072B2",
    "orange": "#D55E00",
    "green": "#009E73",
    "purple": "#6A3D9A",
    "grey": "#6B7280",
    "light": "#E5E7EB",
}

mpl.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "axes.titlesize": 11,
        "axes.labelsize": 9,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "figure.dpi": 160,
        "savefig.dpi": 320,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)


def panel_label(ax, label):
    ax.text(
        -0.13,
        1.08,
        label,
        transform=ax.transAxes,
        fontsize=12,
        fontweight="bold",
        va="top",
    )


def save(fig, stem):
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.png", bbox_inches="tight")
    plt.close(fig)


def correlation_landscape(input_dir):
    data = pd.read_csv(input_dir / "methylation_expression_correlation_data.csv")
    data["adj_p_value"] = pd.to_numeric(data["adj_p_value"], errors="coerce")
    data["correlation"] = pd.to_numeric(data["correlation"], errors="coerce")
    data["significant"] = data["adj_p_value"] < 0.05
    data["neglog10_fdr"] = -np.log10(data["adj_p_value"].clip(lower=1e-300))

    sig = data.loc[data["significant"]].copy()
    region_order = ["upstream", "gene_body", "downstream"]
    region_labels = ["Upstream", "Gene body", "Downstream"]
    counts = (
        sig.groupby("region", observed=True)
        .size()
        .reindex(region_order, fill_value=0)
    )

    fig, axes = plt.subplots(1, 3, figsize=(10.4, 3.1), gridspec_kw={"width_ratios": [2.0, 1.0, 1.25]})

    nonsig = data.loc[~data["significant"]]
    axes[0].scatter(
        nonsig["correlation"],
        nonsig["neglog10_fdr"],
        s=8,
        color=COLORS["light"],
        linewidth=0,
        alpha=0.55,
        rasterized=True,
        label="FDR ≥ 0.05",
    )
    axes[0].scatter(
        sig["correlation"],
        sig["neglog10_fdr"],
        s=15,
        c=np.where(sig["correlation"] >= 0, COLORS["blue"], COLORS["orange"]),
        linewidth=0,
        alpha=0.85,
        rasterized=True,
        label="FDR < 0.05",
    )
    axes[0].axhline(-np.log10(0.05), color=COLORS["grey"], linestyle="--", linewidth=0.8)
    axes[0].axvline(0, color=COLORS["grey"], linewidth=0.7)
    axes[0].set(xlabel="Methylation–isoform correlation", ylabel=r"$-\log_{10}$(FDR)")
    axes[0].set_title("Archived association landscape")
    panel_label(axes[0], "a")

    axes[1].bar(
        range(3),
        counts.values,
        color=[COLORS["grey"], COLORS["blue"], COLORS["grey"]],
        width=0.68,
    )
    axes[1].set_xticks(range(3), region_labels, rotation=25, ha="right")
    axes[1].set_ylabel("FDR-significant rows")
    axes[1].set_title("Genomic localization")
    for i, value in enumerate(counts.values):
        axes[1].text(i, value + max(counts.max() * 0.025, 1), str(value), ha="center", fontsize=8)
    panel_label(axes[1], "b")

    top = (
        sig.assign(abs_r=sig["correlation"].abs())
        .sort_values(["adj_p_value", "abs_r"], ascending=[True, False])
        .drop_duplicates("gene_id")
        .head(12)
        .sort_values("correlation")
    )
    y = np.arange(len(top))
    colors = np.where(top["correlation"] >= 0, COLORS["blue"], COLORS["orange"])
    axes[2].hlines(y, 0, top["correlation"], color=colors, linewidth=1.4)
    axes[2].scatter(top["correlation"], y, c=colors, s=28, zorder=3)
    axes[2].axvline(0, color=COLORS["grey"], linewidth=0.7)
    axes[2].set_yticks(y, top["gene_id"].str.capitalize())
    axes[2].set_xlim(-1.05, 1.05)
    axes[2].set_xlabel("Correlation")
    axes[2].set_title("Top unique genes")
    panel_label(axes[2], "c")

    fig.suptitle(
        "Forebrain-only associations are concentrated in gene bodies",
        fontsize=12,
        fontweight="bold",
        y=1.04,
    )
    fig.tight_layout()
    save(fig, "figure1_correlation_landscape")


def sequence_chromatin_context(input_dir):
    kmers = pd.read_csv(input_dir / "significant_enriched_kmers_categorized.csv")
    kmers["Log2Enr"] = pd.to_numeric(kmers["Log2Enr"], errors="coerce")
    kmers["NegLog10Padj"] = pd.to_numeric(kmers["NegLog10Padj"], errors="coerce")
    kmers = kmers.sort_values(["NegLog10Padj", "Log2Enr"], ascending=False).head(12)
    kmers = kmers.sort_values("Log2Enr")

    hmm = pd.read_csv(input_dir / "enrichment_analysis_results.csv")
    hmm["adj_p_value"] = pd.to_numeric(hmm["adj_p_value"], errors="coerce")
    hmm["log2_fold_change"] = pd.to_numeric(hmm["log2_fold_change"], errors="coerce")
    hmm_sig = hmm.loc[hmm["adj_p_value"] < 0.05].copy()
    hmm_sig["score"] = hmm_sig["log2_fold_change"].abs()
    hmm_top = (
        hmm_sig.sort_values("score", ascending=False)
        .drop_duplicates(["switch_type", "hmm_state"])
        .head(16)
        .sort_values("log2_fold_change")
    )

    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.2))

    y = np.arange(len(kmers))
    kcolors = [
        COLORS["purple"] if "CG" in k else COLORS["grey"] for k in kmers["Kmer"]
    ]
    axes[0].barh(y, kmers["Log2Enr"], color=kcolors, alpha=0.9)
    axes[0].set_yticks(y, kmers["Kmer"])
    axes[0].set_xlabel(r"$\log_2$ enrichment")
    axes[0].set_title("Leading enriched 7-mers")
    axes[0].text(
        0.98,
        0.03,
        f"{kmers['Kmer'].str.contains('CG').sum()}/{len(kmers)} contain CG",
        transform=axes[0].transAxes,
        ha="right",
        fontsize=8,
        color=COLORS["purple"],
    )
    panel_label(axes[0], "a")

    labels = hmm_top["switch_type"].astype(str) + " · " + hmm_top["hmm_state"].astype(str)
    y2 = np.arange(len(hmm_top))
    hcolors = np.where(hmm_top["log2_fold_change"] >= 0, COLORS["green"], COLORS["orange"])
    axes[1].barh(y2, hmm_top["log2_fold_change"], color=hcolors, alpha=0.9)
    axes[1].set_yticks(y2, labels)
    axes[1].axvline(0, color=COLORS["grey"], linewidth=0.7)
    axes[1].set_xlabel(r"$\log_2$ fold enrichment")
    axes[1].set_title("Largest FDR-significant chromatin contrasts")
    panel_label(axes[1], "b")

    fig.suptitle(
        "Sequence and chromatin context nominate loci but do not establish mechanism",
        fontsize=12,
        fontweight="bold",
        y=1.02,
    )
    fig.tight_layout()
    save(fig, "figure2_sequence_chromatin")


def shape_counteraudit(input_dir):
    shape = pd.read_csv(input_dir / "combined_shape_test_results.csv")
    shape["wilcox_p_value"] = pd.to_numeric(shape["wilcox_p_value"], errors="coerce")
    unique = (
        shape.drop_duplicates(["tissue", "region", "shape_type"])
        .loc[:, ["tissue", "region", "shape_type", "wilcox_p_value"]]
        .copy()
    )
    unique = unique.sort_values("wilcox_p_value")
    n = len(unique)
    ranked = np.arange(1, n + 1)
    raw_bh = unique["wilcox_p_value"].to_numpy() * n / ranked
    unique["bh_fdr"] = np.minimum.accumulate(raw_bh[::-1])[::-1]
    unique["bh_fdr"] = unique["bh_fdr"].clip(upper=1)
    unique["label"] = (
        unique["tissue"].str.replace("brain", "", case=False)
        + " · "
        + unique["region"]
        + " · "
        + unique["shape_type"]
    )

    fig, ax = plt.subplots(figsize=(8.0, 4.3))
    y = np.arange(n)
    ax.scatter(unique["bh_fdr"], y, color=COLORS["grey"], s=28)
    ax.axvline(0.05, color=COLORS["orange"], linestyle="--", linewidth=1)
    ax.set_yticks(y, unique["label"])
    ax.set_xlim(0, 1.02)
    ax.invert_yaxis()
    ax.set_xlabel("BH-adjusted Wilcoxon p-value across 16 unique tests")
    ax.set_title(
        "DNA-shape comparisons do not survive family-wise multiplicity control",
        fontweight="bold",
    )
    ax.text(
        0.98,
        0.03,
        "0/16 at FDR < 0.05",
        transform=ax.transAxes,
        ha="right",
        fontsize=9,
        color=COLORS["orange"],
    )
    fig.tight_layout()
    save(fig, "figure3_shape_counteraudit")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Regenerate the archive-level figures.")
    parser.add_argument(
        "--input-dir",
        type=Path,
        help=(
            "Directory containing the four archived CSV inputs. Defaults to "
            "METHYLATION_INPUT_DIR or the historical workspace root."
        ),
    )
    arguments = parser.parse_args()
    input_dir = (arguments.input_dir or DEFAULT_INPUT_DIR).resolve()
    correlation_landscape(input_dir)
    sequence_chromatin_context(input_dir)
    shape_counteraudit(input_dir)
