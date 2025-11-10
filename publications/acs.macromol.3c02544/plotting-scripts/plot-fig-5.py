#!/usr/bin/env python


import os

import numpy as np
import pandas as pd
from pylimer_doctorate_utils.matplotlibUtils import equalify_axes_limits
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

col_mmt_extracted = "G_MMT_MD [MPa]"

data_means = data.groupby(
    ["M", "r sharaf", "w_s sharaf", "row sharaf"], as_index=False
).mean(numeric_only=True)

data_per_col = data_means.groupby(["M"], as_index=False)

# First subplot (a)
for m, group in data_per_col:
    m_clean = group["M"].iloc[0]
    util.scatter_with_mass_style(
        ax1,
        group["G_MD [MPa]"],
        group[col_mmt_extracted],
        m_clean,
        zorder=10,
    )

ax1.set_xlabel("$G_{{\\mathrm{{eq, MD}}}}$ [MPa]")
ax1.set_ylabel("$G_{{\\mathrm{{eq, MD-MMT}}}}$ [MPa]")
ax1.set(xscale="log", yscale="log")
# set limit to max at 1
ax1.set(
    xlim=[0.03, 0.4],
    ylim=[0.03, 0.4],
)
# set tick numbers at 0.03, 0.1 and 0.4
tick_values = [0.03, 0.1, 0.4]
tick_labels = ["0.03", "0.1", "0.4"]
ax1.set_xticks(tick_values)
ax1.set_xticklabels(tick_labels)
ax1.set_yticks(tick_values)
ax1.set_yticklabels(tick_labels)

# equalify_axes_limits(ax)
ax1.plot(ax1.get_xlim(), ax1.get_ylim(), color="black", linestyle="dotted", linewidth=1)
# Second subplot (b)
for m, group in data_per_col:
    m_clean = group["M"].iloc[0]
    util.scatter_with_mass_style(
        ax2,
        group["G_MD [MPa]"],
        group[col_mmt_extracted] / group["G_MD [MPa]"],
        m_clean,
        zorder=10,
    )

ratio = data_means[col_mmt_extracted] / data_means["G_MD [MPa]"]
ax2.set(
    xlim=[0, 0.4],
    ylim=[0, 3.0],
)

ax2.axhline(np.mean(ratio), color="black", linestyle="dashed", linewidth=1.0)  # type: ignore
# shade region around np.mean(ratio) with standard deviation of ratio
std_dev = np.std(ratio)
ax2.fill_between(
    ax2.get_xlim(),
    np.mean(ratio) - std_dev,
    np.mean(ratio) + std_dev,
    color="lightgray",
)

ax2.set_xlabel("$G_{{\\mathrm{{eq, MD}}}}$ [MPa]")
ax2.set_ylabel("$G_{{\\mathrm{{eq, MD-MMT}}}} / G_{{\\mathrm{{eq, MD}}}}$ [MPa]")

# Add subplot labels and shared legend
util.add_subplot_labels(axes)
util.setup_shared_legend(fig, axes)


saveFigure(fig, os.path.join(base_dir, "plots", "plot-fig-5"))
dismissFigure(fig)
