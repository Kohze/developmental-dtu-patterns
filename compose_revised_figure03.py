#!/usr/bin/env python3
"""Compose main Figure 3 from existing calibrated and supplementary assets."""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
from tempfile import TemporaryDirectory

import matplotlib as mpl
import matplotlib.pyplot as plt
from pypdf import PdfReader, PdfWriter, Transformation
from pypdf._page import PageObject


ROOT = Path(__file__).resolve().parent
DTU = ROOT / "dtu_analysis" / "figures"
TOP_SOURCE = DTU / "figure6_transient_regional_episodes.pdf"
BOTTOM_SOURCE = ROOT / "figures" / "figureS6_repaired_candidate_trajectories.pdf"
OUTPUT = ROOT / "figures" / "figure03_repaired_transient_episodes.pdf"

PAGE_WIDTH = 518.4
TOP_HEIGHT = 198.0
PANEL_GAP = 21.0


def load_page(path: Path):
    return PdfReader(str(path)).pages[0]


top = load_page(TOP_SOURCE)
bottom = load_page(BOTTOM_SOURCE)

top_width = float(top.mediabox.width)
top_full_height = float(top.mediabox.height)
bottom_width = float(bottom.mediabox.width)
bottom_height = float(bottom.mediabox.height)
bottom_scale = PAGE_WIDTH / bottom_width
scaled_bottom_height = bottom_height * bottom_scale
page_height = TOP_HEIGHT + PANEL_GAP + scaled_bottom_height

canvas = PageObject.create_blank_page(width=PAGE_WIDTH, height=page_height)

# Retain panel A from the calibrated scan. PDF coordinates start at the lower
# left, so the source strip is translated down into the new page.
top.mediabox.lower_left = (0, top_full_height - TOP_HEIGHT)
top.mediabox.upper_right = (top_width, top_full_height)
top.cropbox.lower_left = top.mediabox.lower_left
top.cropbox.upper_right = top.mediabox.upper_right
canvas.merge_transformed_page(
    top,
    Transformation().translate(
        tx=0,
        ty=scaled_bottom_height + PANEL_GAP - (top_full_height - TOP_HEIGHT),
    ),
)

canvas.merge_transformed_page(
    bottom,
    Transformation().scale(bottom_scale).translate(tx=0, ty=0),
)

# Add the panel label without rasterising either source graph.
with TemporaryDirectory(prefix="figure03-label-") as temporary:
    overlay_path = Path(temporary) / "label.pdf"
    mpl.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42})
    fig = plt.figure(figsize=(PAGE_WIDTH / 72, page_height / 72))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.axis("off")
    ax.text(
        0.013,
        scaled_bottom_height / page_height - 0.006,
        "B",
        fontsize=11,
        fontweight="bold",
        va="top",
        color="#25313C",
    )
    fig.savefig(
        overlay_path,
        transparent=True,
        metadata={
            "CreationDate": datetime(2000, 1, 1, tzinfo=timezone.utc),
            "ModDate": datetime(2000, 1, 1, tzinfo=timezone.utc),
        },
    )
    plt.close(fig)
    canvas.merge_page(PdfReader(str(overlay_path)).pages[0])

writer = PdfWriter()
writer.add_page(canvas)
writer.add_metadata(
    {
        "/Title": "Calibrated transient episodes and reconstructed legacy trajectories",
        "/CreationDate": "D:20000101000000Z",
        "/ModDate": "D:20000101000000Z",
    }
)
with OUTPUT.open("wb") as handle:
    writer.write(handle)

print(OUTPUT)
