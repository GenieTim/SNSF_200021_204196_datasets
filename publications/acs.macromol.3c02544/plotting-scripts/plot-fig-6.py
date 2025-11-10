#!/usr/bin/env python
import os

import pandas as pd
from pylimer_doctorate_utils.saveFigure import dismissFigure, saveFigure

import util

base_dir = os.path.join(os.path.dirname(__file__))

data = pd.read_csv(
    os.path.join(base_dir, "all-structures-r-sharaf-wsol-sharaf-without-wsol.csv")
)
data = data[(data["M"] == 7380) | (data["M"] == 7500) | (data["M"] == 10800)]

# Setup subplot figure
fig, axes = util.setup_ab_figure()
ax1, ax2 = axes["A"], axes["B"]

data_means = data.groupby(
    ["M", "r sharaf", "w_s sharaf", "row sharaf"], as_index=False
).mean(numeric_only=True)

data_per_col = data_means.groupby(["M"], as_index=False)

# First subplot (a)
for m, group in data_per_col:
    m_clean = group["M"].iloc[0]
    util.scatter_with_mass_style(
        ax1,
        group["G_ANT [MPa]"],
        group["G_MMT_MD_phantom [MPa]"],
        m_clean,
        zorder=10,
    )

ax1.set_xlabel("$G_{{\\mathrm{{eq, ANT}}}}$ [MPa]")
ax1.set_ylabel("$G_{{\\mathrm{{eq, MD-MMT,junctions}}}}$ [MPa]")
ax1.set(xscale="log", yscale="log")
# set limit to max at 1
ax1.set(
    xlim=[0.03, 0.25],
    ylim=[0.03, 0.25],
)
# set tick numbers at 0.03, 0.1 and 0.25
tick_values = [0.03, 0.1, 0.25]
ax1.set_xticks(tick_values)
ax1.set_yticks(tick_values)
ax1.tick_params(which="minor", labelbottom=False, labelleft=False)
# ax.ticklabel_format(style="plain", axis="both", useOffset=False)
# ax.ticklabel_format(scilimits=(-3, 3))
# set tick formatter to scalar
ax1.xaxis.set_major_formatter("{x:.2f}")
ax1.yaxis.set_major_formatter("{x:.2f}")

# equalify_axes_limits(ax)
ax1.plot(ax1.get_xlim(), ax1.get_ylim(), color="black", linestyle="dotted", linewidth=1)
# Second subplot (b)
for m, group in data_per_col:
    m_clean = group["M"].iloc[0]
    util.scatter_with_mass_style(
        ax2,
        group["G_MD [MPa]"] - group["G_MMT_MD_phantom [MPa]"],
        group["G_MMT_MD [MPa]"] - group["G_MMT_MD_phantom [MPa]"],
        m_clean,
        zorder=10,
    )

ax2.set(xscale="log", yscale="log")
# set limit to max at 1
ax2.set(
    xlim=[0.03, 0.25],
    ylim=[0.03, 0.25],
)
# set tick numbers at 0.03, 0.1 and 0.25
tick_values = [0.03, 0.1, 0.25]
ax2.set_xticks(tick_values)
ax2.set_yticks(tick_values)
ax2.tick_params(which="minor", labelbottom=False, labelleft=False)
# ax.ticklabel_format(style="plain", axis="both", useOffset=False)
ax2.xaxis.set_major_formatter("{x:.2f}")
ax2.yaxis.set_major_formatter("{x:.2f}")

# equalify_axes_limits(ax)
ax2.plot(ax2.get_xlim(), ax2.get_ylim(), color="black", linestyle="dotted", linewidth=1)

ax2.set_xlabel("$G_{{\\mathrm{{eq, MD,entanglements}}}}$ [MPa]")
ax2.set_ylabel("$G_{{\\mathrm{{eq, MD-MMT,entanglements}}}}$ [MPa]")

# Add subplot labels and shared legend
util.add_subplot_labels(axes)
util.setup_shared_legend(fig, axes)

saveFigure(fig, os.path.join(base_dir, "plots", "plot-fig-6"))
dismissFigure(fig)
