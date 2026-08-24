"""Create the candidate architecture and complementary-evidence figure.

Panel A is a RefSeq interpretation layer. Panel B compares independent
prioritisation axes within the archived analysis. Its factual inputs are
versioned in tables/candidate_mechanism_crosswalk.csv and
tables/candidate_counter_audit.csv.
"""

from pathlib import Path
import csv
from datetime import datetime, timezone

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "tables" / "candidate_mechanism_crosswalk.csv"
COUNTER_SOURCE = ROOT / "tables" / "candidate_counter_audit.csv"
OUT_PDF = ROOT / "figures" / "figure7_candidate_mechanism.pdf"
OUT_PNG = ROOT / "figures" / "figure7_candidate_mechanism.png"
PDF_METADATA = {
    "CreationDate": datetime(2000, 1, 1, tzinfo=timezone.utc),
    "ModDate": datetime(2000, 1, 1, tzinfo=timezone.utc),
}

mpl.rcParams.update(
    {
        "font.family": "Arial",
        "font.size": 8.8,
        "axes.linewidth": 0.8,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "savefig.facecolor": "white",
    }
)

with SOURCE.open(encoding="utf-8", newline="") as handle:
    rows = list(csv.DictReader(handle))
with COUNTER_SOURCE.open(encoding="utf-8", newline="") as handle:
    counter_rows = list(csv.DictReader(handle))

by_gene = {row["gene"]: row for row in rows}

BLUE = "#0072B2"
ORANGE = "#D55E00"
GREEN = "#009E73"
PURPLE = "#CC79A7"
DARK = "#25313C"
MID = "#60717F"
LIGHT = "#F2F5F7"
PALE_BLUE = "#E8F3F8"
PALE_ORANGE = "#FCEFE9"
PALE_GREEN = "#E8F5F0"
PALE_PURPLE = "#F7EAF2"

fig = plt.figure(figsize=(7.194, 6.35))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.axis("off")

ax.text(0.035, 0.965, "A", fontsize=11, fontweight="bold", color=DARK, va="top")
ax.text(
    0.075,
    0.965,
    "RefSeq transcript architecture and direction of the E15.5 exchange",
    fontsize=9.2,
    fontweight="bold",
    color=DARK,
    va="top",
)

architecture = [
    ("Ntrk2", "truncated TrkB.T1", "full-length TrkB.FL + kinase", BLUE),
    ("Gpm6a", "short, distinct N-terminus", "long isoform", ORANGE),
    ("Bin1", "short central architecture", "588-aa isoform", ORANGE),
    ("Tecr", "long isoform", "in-frame exon absent", GREEN),
    ("Scg3", "long isoform", "internal in-frame shortening", PURPLE),
    ("Armc8", "long + HEAT annotation", "short, distinct C-terminus", MID),
]

start_y = 0.895
row_gap = 0.095
for index, (gene, higher_label, lower_label, colour) in enumerate(architecture):
    y = start_y - index * row_gap
    ax.text(0.055, y, gene, fontstyle="italic", fontweight="bold", color=DARK, va="center")

    low_box = FancyBboxPatch(
        (0.185, y - 0.025),
        0.265,
        0.05,
        boxstyle="round,pad=0.006,rounding_size=0.012",
        facecolor=LIGHT,
        edgecolor="#A9B4BC",
        linewidth=0.9,
    )
    high_box = FancyBboxPatch(
        (0.57, y - 0.025),
        0.35,
        0.05,
        boxstyle="round,pad=0.006,rounding_size=0.012",
        facecolor="white",
        edgecolor=colour,
        linewidth=1.5,
    )
    ax.add_patch(low_box)
    ax.add_patch(high_box)

    row = by_gene[gene]
    ax.text(
        0.198,
        y + 0.008,
        row["lower_midbrain_accession"],
        fontsize=7.5,
        fontweight="bold",
        color=DARK,
        va="center",
    )
    ax.text(0.198, y - 0.011, lower_label, fontsize=6.8, color=MID, va="center")
    ax.text(
        0.585,
        y + 0.008,
        row["higher_midbrain_accession"],
        fontsize=7.5,
        fontweight="bold",
        color=colour,
        va="center",
    )
    ax.text(0.585, y - 0.011, higher_label, fontsize=6.8, color=DARK, va="center")

    arrow = FancyArrowPatch(
        (0.465, y),
        (0.552, y),
        arrowstyle="-|>",
        mutation_scale=11,
        linewidth=1.3,
        color=colour,
    )
    ax.add_patch(arrow)
    ax.text(0.508, y + 0.021, "E15.5", fontsize=6.2, color=MID, ha="center")

ax.plot([0.035, 0.965], [0.305, 0.305], color="#D9E0E5", linewidth=1)

ax.text(0.035, 0.275, "B", fontsize=11, fontweight="bold", color=DARK, va="top")
ax.text(
    0.075,
    0.275,
    "Complementary evidence axes prioritise distinct candidates",
    fontsize=9.2,
    fontweight="bold",
    color=DARK,
    va="top",
)

scan_order = sorted(counter_rows, key=lambda row: int(row["calibrated_scan_rank"]))
effect_order = sorted(counter_rows, key=lambda row: int(row["effect_size_rank"]))
joint_order = sorted(counter_rows, key=lambda row: int(row["joint_model_rank_of_4577"]))

audit_boxes = [
    (
        0.055,
        0.158,
        0.425,
        0.072,
        PALE_PURPLE,
        PURPLE,
        "PATTERN-DETECTOR RANK",
        f"{scan_order[0]['gene']}  #1   |   {scan_order[1]['gene']}  #2",
        "best globally adjusted diverge--reconverge evidence",
    ),
    (
        0.52,
        0.158,
        0.425,
        0.072,
        PALE_ORANGE,
        ORANGE,
        "EFFECT SIZE",
        (
            f"{effect_order[0]['gene']}  "
            f"{float(effect_order[0]['max_abs_usage_difference']):.3f}"
            f"   |   {effect_order[1]['gene']}  "
            f"{float(effect_order[1]['max_abs_usage_difference']):.3f}"
        ),
        "largest midbrain-versus-comparison fraction differences",
    ),
    (
        0.055,
        0.060,
        0.425,
        0.072,
        PALE_GREEN,
        GREEN,
        "FOCUSED COUNT MODEL",
        (
            f"{joint_order[0]['gene']}  #{joint_order[0]['joint_model_rank_of_4577']}"
            f"   |   {joint_order[1]['gene']}  #{joint_order[1]['joint_model_rank_of_4577']}"
        ),
        "Ntrk2 also passes gene q=0.024 in the shared-sample model",
    ),
    (
        0.52,
        0.060,
        0.425,
        0.072,
        PALE_BLUE,
        BLUE,
        "PROTEIN INTERPRETABILITY",
        "Ntrk2: kinase-containing FL versus truncated T1",
        "clearest domain contrast; effect-size rank 6 of 6",
    ),
]

for x, y, w, h, face, edge, title, finding, qualifier in audit_boxes:
    box = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.006,rounding_size=0.014",
        facecolor=face,
        edgecolor=edge,
        linewidth=1.2,
    )
    ax.add_patch(box)
    ax.text(x + 0.015, y + h - 0.019, title, fontsize=6.0, fontweight="bold", color=edge)
    ax.text(x + 0.015, y + 0.031, finding, fontsize=6.9, fontweight="bold", color=DARK)
    ax.text(x + 0.015, y + 0.011, qualifier, fontsize=5.6, color=MID)

ax.text(
    0.5,
    0.021,
    "The six-gene panel spans transcript architecture, effect size and model support.",
    fontsize=7.0,
    fontweight="bold",
    color=DARK,
    ha="center",
)

OUT_PDF.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT_PDF, bbox_inches=None, metadata=PDF_METADATA)
fig.savefig(OUT_PNG, dpi=300, bbox_inches=None)
plt.close(fig)
