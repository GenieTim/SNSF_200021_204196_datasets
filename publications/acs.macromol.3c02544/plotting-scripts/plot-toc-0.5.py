#!/usr/bin/env python
import os

import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.ticker import MultipleLocator
from pylimer_doctorate_utils.saveFigure import addMinorTicks, saveFigure, setupAxis


base_dir = os.path.join(os.path.dirname(__file__))

data = pd.read_excel(os.path.join(base_dir, "manual_table.xlsx"))
data.sort_values(by="G_MD [MPa]", inplace=True)

# fig, ax = setupFigure(8, 4.5) if not os.environ.get("FIGURE_MODE", "") == "thesis" else setupFigure()
fig, ax_dict = plt.subplot_mosaic(
    "L;A", figsize=(7, 7), layout="constrained", height_ratios=[1, 7]
)
fig.get_layout_engine().set(w_pad=0, h_pad=0, hspace=0.0175, wspace=0.0175)
ax_dict["L"].axis("off")  # Hide the legend axis

ax = ax_dict["A"]
setupAxis(ax)

ax2 = ax.twinx()
addMinorTicks(ax2)

p1 = ax.scatter(
    data["G_MD [MPa]"],
    data["G_MMT_MD [MPa]"] / data["G_MD [MPa]"],
    s=60,
    label="MD: this work",
)
# ax2.scatter([], [])  # skip color
p2 = ax.scatter(
    data["G_MD [MPa]"],
    data["G_EXP [MPa]"] / data["G_MD [MPa]"],
    marker="D",
    s=60,
    label="Exp.",
)
ax.axhline(1.0, color="black", linewidth=1.0, linestyle="dashed")
xlim = ax.get_xlim()
# ax2.fill_between([0., 1.], [
#                  0.9, 0.9], [1.1, 1.1], alpha=0.25)

# labelinging
ax.set_xlabel("$G_{{\\mathrm{{eq, MD}}}}$ [MPa]")
ax.set_ylabel(
    "$G_{{\\mathrm{{eq, MMT}}}}/G_{{\\mathrm{{eq, MD}}}}$", color=p1.get_facecolor()
)
ax2.set_ylabel(
    "$G_{{\\mathrm{{eq, Exp}}}}/G_{{\\mathrm{{eq, MD}}}}$", color=p2.get_facecolor()
)  # , rotation=270, labelpad=20)

# ax.text(
#     0.015,
#     0.325,
#     "Exp.: Sharaf et al., 1996\nMD & MMT: this work",
#     verticalalignment="top",
#     bbox={"boxstyle": "square", "alpha": 0.5, "facecolor": "white"},
#     linespacing=1.75,
# )

plt.xticks([0.1 * i for i in range(10)])

xlim = [0.0, 0.4]  # xlim[1]]
ylim = [0.0, 1.401]
ax.set(ylim=ylim, xlim=xlim)
ax2.set(ylim=ylim, xlim=xlim)

ax.get_xaxis().set_minor_locator(MultipleLocator(0.05))
ax.get_yaxis().set_minor_locator(MultipleLocator(0.1))
ax2.get_yaxis().set_minor_locator(MultipleLocator(0.1))


# for item in ax.get_xticklabels() + ax.get_yticklabels() + ax2.get_yticklabels():
#     item.set_fontsize(16)
# for item in [ax.xaxis.label, ax.yaxis.label, ax2.yaxis.label]:
#     item.set_fontsize(20)

# ax.legend(loc="lower left")

# Collect all handles and labels from axes A
handles, labels = [], []
for handle, label in zip(*ax_dict["A"].get_legend_handles_labels()):
    handles.append(handle)
    labels.append(label)

# Create an empty handle for spacing
empty_handle = Rectangle((0, 0), 1, 1, fill=False, edgecolor="none", visible=False)

centered_handles = [
    *[empty_handle for _ in range(len(handles))],
    *handles,
    *[empty_handle for _ in range(len(handles))],
]
centered_labels = [
    *["" for _ in range(len(labels))],
    *labels,
    *["" for _ in range(len(handles))],
]

lgd = ax_dict["L"].legend(
    # handles,
    # labels,
    centered_handles,
    centered_labels,
    loc="center",
    mode="expand",
    ncol=3,
    frameon=True,
    fontsize=8,
    borderaxespad=0.0,
    columnspacing=0.0,  # Reduce space between columns (default is 2.0)
    handletextpad=0.4,
    # Reduce space between symbol and text (default is 0.8)
)
# add space on top of the figure for the new legend
fig.set_figheight(fig.get_figheight() + lgd.get_window_extent().height / fig.dpi + 0.45)

saveFigure(fig, os.path.join(base_dir, "plots", "toc-2sided-v2"), ignore_legend=True)
