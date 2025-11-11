#!/usr/bin/env python


import os

import numpy as np
import pandas as pd
from pint import UnitRegistry
from pylimer_doctorate_utils.saveFigure import dismissFigure, saveFigure
from pylimer_tools.calc.miller_macosko_theory import (predict_p_from_w_sol,
                                                      predict_shear_modulus)

import util

base_dir = os.path.join(os.path.dirname(__file__))

data = pd.read_csv(
    os.path.join(base_dir, "..", "processed-data", "all-structures-r-sharaf-wsol-sharaf-without-wsol.csv")
)

# Setup subplot figure
fig, axes = util.setup_ab_figure()
ax1, ax2 = axes["A"], axes["B"]

ureg = UnitRegistry()

# where the column "r sharaf" is missing,
# replace it with the values from "r_x" or "r_y"
data["r sharaf"] = data["r sharaf"].fillna(data["r_x"].fillna(data["r_y"]))

col_mmt_extracted = "G_MMT [MPa] 2"

data_means = data.groupby(
    ["M", "r sharaf", "w_s sharaf", "row sharaf"], as_index=False
).mean(numeric_only=True)

del data

data_means[col_mmt_extracted] = data_means.apply(
    lambda row: predict_shear_modulus(
        ureg=ureg,
        r=row["r sharaf"],
        p=predict_p_from_w_sol(
            w_sol=row["w_s sharaf"],
            r=row["r sharaf"],
            functionality_per_type={1: 2, 2: 4},
            weight_fractions={1: 1, 2: 0},
            b2=1.0,
        ),
        functionality_per_type={1: 2, 2: 4},
        crosslinker_type=2,
        f=4,
        temperature=300.0 * ureg("K"),
        g_e_1=0.207 * ureg("MPa"),
        b2=1.0,
        nu=row["nu [nm^-3]"] * ureg("nm^-3"),
        network=None,
    )
    .to("MPa")
    .magnitude,
    axis=1,
)

data_per_col = data_means.groupby(["M"], as_index=False)

# First subplot (a)
for m, group in data_per_col:
    m_clean = group["M"].iloc[0]
    util.scatter_with_mass_style(
        ax1,
        group["G [N mm^-2] sharaf"],
        group[col_mmt_extracted],
        m_clean,
        zorder=10,
    )

ax1.set_xlabel("$G_{{\\mathrm{{eq, exp, extracted}}}}$ [MPa]")
ax1.set_ylabel("$G_{{\\mathrm{{eq, MMT, extracted}}}}$ [MPa]")
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

ax1.plot(ax1.get_xlim(), ax1.get_ylim(), color="black", linestyle="dotted", linewidth=1)
# Second subplot (b)
for m, group in data_per_col:
    m_clean = group["M"].iloc[0]
    util.scatter_with_mass_style(
        ax2,
        group["G [N mm^-2] sharaf"],
        group[col_mmt_extracted] / group["G [N mm^-2] sharaf"],
        m_clean,
        zorder=10,
    )

ax2.set_ylabel(
    "$G_{{\\mathrm{{eq, MMT, extracted}}}} / G_{{\\mathrm{{eq, exp, extracted}}}}$ [MPa]"
)
ax2.set_xlabel("$G_{{\\mathrm{{eq, exp, extracted}}}}$ [MPa]")

ax2.set(
    xlim=[0, 0.3],
    ylim=[0, 3.0],
)

ratio = data_means[col_mmt_extracted] / data_means["G [N mm^-2] sharaf"]
# drop the highest value
ratio = ratio[ratio < 2.0]

ax2.axhline(np.mean(ratio), color="black", linestyle="dashed", linewidth=1.0)  # type: ignore
# shade region around np.mean(ratio) with standard deviation of ratio
std_dev = np.std(ratio)
ax2.fill_between(
    ax2.get_xlim(),
    np.mean(ratio) - std_dev,
    np.mean(ratio) + std_dev,
    color="lightgray",
)

# Add subplot labels and shared legend
util.add_subplot_labels(axes)
util.setup_shared_legend(fig, axes, ncol=5)


saveFigure(fig, os.path.join(base_dir, "plots", "plot-fig-2"))
dismissFigure(fig)
