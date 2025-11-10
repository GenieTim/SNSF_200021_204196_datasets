#!/usr/bin/env python


import glob
import math
import os
import pickle
import random
import warnings

import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
from pint import UnitRegistry
from pylimer_doctorate_utils.fileConfigResolver import (
    extractDataFromFilename, findStructureFileForOutput,
    fullyAnalyseStructureFile)
from pylimer_doctorate_utils.matplotlibUtils import equalify_axes_limits
from pylimer_doctorate_utils.progressbar import progressbar
from pylimer_doctorate_utils.saveFigure import (addAxisClone, addMinorTicks,
                                                correctFontSizes,
                                                default_font_size,
                                                dismissFigure, saveFigure,
                                                setupFigure)
from pylimer_tools.calc.miller_macosko_theory import (
    compute_modulus_decomposition, predict_p_from_w_sol)

import util

base_dir = os.path.join(os.path.dirname(__file__))

data = pd.read_csv(
    os.path.join(base_dir, "all-structures-r-sharaf-wsol-sharaf-without-wsol.csv")
)

fig, ax = setupFigure()

ureg = UnitRegistry()

data_means = data.groupby(
    ["M", "r sharaf", "w_s sharaf", "row sharaf"], as_index=False
).mean(numeric_only=True)

data_per_col = data_means.groupby(["M"], as_index=False)
for m, group in data_per_col:
    m_clean = group["M"].iloc[0]

    def compute_xy_values(row):
        g_mmt_phantom, g_mmt_entanglement, _, _ = compute_modulus_decomposition(
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
        )

        return (
            g_mmt_phantom.to("MPa").magnitude,
            g_mmt_entanglement.to("MPa").magnitude,
        )

    xy_values = [compute_xy_values(row) for _, row in group.iterrows()]
    y_values = [xy[0] for xy in xy_values]
    x_values = [xy[1] for xy in xy_values]

    util.scatter_with_mass_style(
        ax,
        x_values,
        y_values,
        m_clean,
    )


ax.set(xscale="log", yscale="log")

ax.set(
    xlim=[0.005, 0.2],
    ylim=[0.005, 0.2],
)
equalify_axes_limits(ax)

# add dotted black diagnonal line
ax.plot(ax.get_xlim(), ax.get_ylim(), color="black", linestyle="dotted", linewidth=1)

# add labels
ax.set_xlabel("$G_{{\\mathrm{{eq, entanglements, extracted}}}}$ [MPa]")
ax.set_ylabel("$G_{{\\mathrm{{eq, junctions, extracted}}}}$ [MPa]")
ax.legend()


saveFigure(fig, os.path.join(base_dir, "plots", "plot-fig-3"))
