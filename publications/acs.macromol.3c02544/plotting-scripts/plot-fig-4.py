#!/usr/bin/env python

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pint import UnitRegistry
from pylimer_doctorate_utils.saveFigure import (dismissFigure, saveFigure,
                                                setupFigure)
from pylimer_tools.calc.miller_macosko_theory import (
    compute_modulus_decomposition, predict_maximum_p)

import util

base_dir = os.path.join(os.path.dirname(__file__))

ureg = UnitRegistry()

fig, ax = setupFigure()

ax.axhline(1.0, color="black", linestyle="--", linewidth=1)

ax.set(xscale="log", yscale="log")
ax.set_xlabel("$M_n$ [g/mol]")
ax.set_ylabel("$G_{{\\mathrm{{eq, entanglement}}}}/G_{{\\mathrm{{eq, junction}}}}$")

# Define colors from matplotlib's tab10 colormap
_tab10_colors = plt.cm.get_cmap("tab10")
color_offset = 5

p = 0.9
for i, r in enumerate([1.3, 1.0, 0.7]):
    mns = np.logspace(2, 6, 10)
    ratios = []
    for mn in mns:
        g_mmt_phantom, g_mmt_entanglement, _, _ = compute_modulus_decomposition(
            ureg=ureg,
            r=r,
            p=p * predict_maximum_p(r=r, f=4, b2=1),
            functionality_per_type={1: 2, 2: 4},
            crosslinker_type=2,
            f=4,
            temperature=300.0 * ureg("K"),
            g_e_1=0.207 * ureg("MPa"),
            b2=1.0,
            nu=0.965
            * ureg("g/cm^3")
            / (mn * ureg("g/mol"))
            * (6.022e23 * ureg("mol^-1")),
        )
        ratios.append(
            g_mmt_entanglement.to("MPa").magnitude / g_mmt_phantom.to("MPa").magnitude
        )
    ax.plot(
        mns, ratios, label=f"$r={r}$", zorder=10, color=_tab10_colors(i + color_offset)
    )

# add an arrow pointing to y = 1., at x = 12000
ax.annotate(
    "",
    xy=(12000, 1.0),
    xytext=(12000, 0.1),
    arrowprops=dict(arrowstyle="->", color="black"),
)

ax.set(
    xlim=(1e2, 1e6),
    ylim=(0.01, 100),
)

ax.legend() 


saveFigure(fig, os.path.join(base_dir, "plots", "plot-fig-4"))
