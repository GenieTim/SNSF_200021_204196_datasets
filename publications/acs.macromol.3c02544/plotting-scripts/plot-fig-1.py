#!/usr/bin/env python

import os

import pandas as pd
from pylimer_doctorate_utils.matplotlibUtils import equalify_axes_limits
from pylimer_doctorate_utils.saveFigure import dismissFigure, saveFigure
from pylimer_tools.calc.miller_macosko_theory import predict_p_from_w_sol

import util

base_dir = os.path.join(os.path.dirname(__file__))

data = pd.read_csv(
    os.path.join(base_dir, "..", "processed-data", "all-structures-r-sharaf-wsol-sharaf-without-wsol.csv")
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
        group["w_s sharaf"],
        group["ws (Fabian)"],
        m_clean,
        zorder=10,
    )

ax1.set_xlabel("$w_{{\\mathrm{{sol, exp}}}}$")
ax1.set_ylabel("$w_{{\\mathrm{{sol, MD}}}}$")
ax1.set(xscale="log", yscale="log")
equalify_axes_limits(ax1)
ax1.set_aspect("equal")
ax1.plot(ax1.get_xlim(), ax1.get_ylim(), color="black", linestyle="dotted", linewidth=1)

# Second subplot (b)
for m, group in data_per_col:
    m_clean = group["M"].iloc[0]
    util.scatter_with_mass_style(
        ax2,
        [
            predict_p_from_w_sol(
                w_sol=row["w_s sharaf"],
                functionality_per_type={
                    1: 2,
                    2: 4,
                },
                r=row["r sharaf"],
                b2=1.0,
                weight_fractions={1: 1, 2: 0},
            )
            for i, row in group.iterrows()
        ],
        group["p"],
        m_clean,
        zorder=10,
    )

ax2.set_xlabel("$p_{{\\mathrm{{exp}}}}$")
ax2.set_ylabel("$p_{{\\mathrm{{MD}}}}$")
equalify_axes_limits(ax2)
ax2.set_aspect("equal")
ax2.plot(ax2.get_xlim(), ax2.get_ylim(), color="black", linestyle="dotted", linewidth=1)

# Add subplot labels and shared legend
util.add_subplot_labels(axes)
util.setup_shared_legend(fig, axes)


saveFigure(fig, os.path.join(base_dir, "plots", "plot-fig-1"))
dismissFigure(fig)
