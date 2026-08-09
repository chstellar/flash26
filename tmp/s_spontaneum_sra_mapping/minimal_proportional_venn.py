#!/usr/bin/env python3
"""
Minimal standalone proportional Venn-style plot.

The large outer circle represents the full set (= 1.0). Two inner circles
represent A and B; their areas are proportional to their shares, and their
overlap area is approximately proportional to overlap_share.
"""

import math
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Circle


# ------------------------- user-editable values -------------------------
out_prefix = Path("proportional_venn")

a_label = "Set A"
b_label = "Set B"
a_share = 0.34
b_share = 0.26
overlap_share = 0.12

title = "Proportional overlap within the full set"
# ------------------------------------------------------------------------


def circle_overlap_area(r1, r2, d):
    if d >= r1 + r2:
        return 0.0
    if d <= abs(r1 - r2):
        return math.pi * min(r1, r2) ** 2
    a = r1**2 * math.acos((d**2 + r1**2 - r2**2) / (2 * d * r1))
    b = r2**2 * math.acos((d**2 + r2**2 - r1**2) / (2 * d * r2))
    c = 0.5 * math.sqrt(
        max(0.0, (-d + r1 + r2) * (d + r1 - r2) * (d - r1 + r2) * (d + r1 + r2))
    )
    return a + b - c


def distance_for_overlap(r1, r2, target_area):
    max_overlap = math.pi * min(r1, r2) ** 2
    target_area = max(0.0, min(target_area, max_overlap))
    lo = abs(r1 - r2)
    hi = r1 + r2
    for _ in range(80):
        mid = (lo + hi) / 2
        if circle_overlap_area(r1, r2, mid) > target_area:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def pct(x):
    return f"{100 * x:.1f}%"


if overlap_share > min(a_share, b_share):
    raise ValueError("overlap_share cannot exceed the smaller of a_share and b_share.")
if a_share + b_share - overlap_share > 1:
    raise ValueError("A union B cannot exceed the full set.")

outer_r = 1.0
r_a = math.sqrt(a_share)
r_b = math.sqrt(b_share)
target_overlap_area = math.pi * overlap_share
d = distance_for_overlap(r_a, r_b, target_overlap_area)

# Keep the two circles centered within the large circle.
x_a = -d / 2
x_b = d / 2
y = 0.0
scale_to_fit = max(
    1.0,
    (abs(x_a) + r_a) / 0.93,
    (abs(x_b) + r_b) / 0.93,
)
r_a /= scale_to_fit
r_b /= scale_to_fit
x_a /= scale_to_fit
x_b /= scale_to_fit

fig, ax = plt.subplots(figsize=(7, 7))
ax.set_aspect("equal")
ax.axis("off")

outer_color = "#F4F1EA"
a_color = "#3B8D7A"
b_color = "#D99A30"
edge_color = "#263238"

ax.add_patch(Circle((0, 0), outer_r, facecolor=outer_color, edgecolor=edge_color, lw=2.2, zorder=1))
ax.add_patch(Circle((x_a, y), r_a, facecolor=a_color, edgecolor=edge_color, lw=1.8, alpha=0.72, zorder=2))
ax.add_patch(Circle((x_b, y), r_b, facecolor=b_color, edgecolor=edge_color, lw=1.8, alpha=0.72, zorder=3))

a_only = a_share - overlap_share
b_only = b_share - overlap_share
neither = 1 - (a_share + b_share - overlap_share)

ax.text(x_a - r_a * 0.30, y + 0.05, f"{a_label}\n{pct(a_only)} only", ha="center", va="center", fontsize=12)
ax.text(x_b + r_b * 0.30, y + 0.05, f"{b_label}\n{pct(b_only)} only", ha="center", va="center", fontsize=12)
ax.text((x_a + x_b) / 2, y - 0.05, f"overlap\n{pct(overlap_share)}", ha="center", va="center", fontsize=12, fontweight="bold")
ax.text(0, -0.88, f"neither: {pct(neither)}", ha="center", va="center", fontsize=12)
ax.text(0, 1.13, title, ha="center", va="center", fontsize=15, fontweight="bold")
ax.text(0, 0.98, "Full set = 100%", ha="center", va="center", fontsize=11, color="#455A64")

ax.set_xlim(-1.18, 1.18)
ax.set_ylim(-1.18, 1.23)

fig.savefig(f"{out_prefix}.png", dpi=300, bbox_inches="tight")
fig.savefig(f"{out_prefix}.pdf", bbox_inches="tight")
print(f"Wrote {out_prefix}.png and {out_prefix}.pdf")
