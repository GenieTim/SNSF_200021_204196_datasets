#!/usr/bin/env python
"""
Utility functions for consistent plotting across all figure scripts.
Provides standardized colors and symbols for different molecular masses.
"""
import os

import matplotlib
import matplotlib.pyplot as plt

# Define the standard masses used in the dataset
STANDARD_MASSES = [7380.0, 7500.0, 10800.0, 15900.0, 23500.0]

# Define colors from matplotlib's tab10 colormap
_tab10_colors = plt.cm.get_cmap("tab10")
MASS_COLORS = {
    7380.0: _tab10_colors(0),  # blue
    7500.0: _tab10_colors(1),  # orange
    10800.0: _tab10_colors(2),  # green
    15900.0: _tab10_colors(3),  # red
    23500.0: _tab10_colors(4),  # purple
}

# Define distinct markers for each mass
MASS_MARKERS = {
    7380.0: "s",  # square
    7500.0: "o",  # circle
    10800.0: "^",  # triangle up
    15900.0: "v",  # triangle down
    23500.0: "D",  # diamond
}


def get_mass_color(mass):
    """
    Get the standardized color for a given molecular mass.

    Parameters
    ----------
    mass : float
        The molecular mass value

    Returns
    -------
    color
        The matplotlib color for the given mass

    Raises
    ------
    ValueError
        If the mass is not in the standard set of masses
    """
    # Round to nearest integer to handle floating point precision
    mass_rounded = round(mass)

    # Find the closest standard mass
    closest_mass = min(STANDARD_MASSES, key=lambda x: abs(x - mass_rounded))

    # Check if it's close enough (within 1% tolerance)
    if abs(closest_mass - mass_rounded) / closest_mass > 0.01:
        raise ValueError(
            f"Mass {mass} is not close to any standard mass: {STANDARD_MASSES}"
        )

    return MASS_COLORS[closest_mass]


def get_mass_marker(mass):
    """
    Get the standardized marker for a given molecular mass.

    Parameters
    ----------
    mass : float
        The molecular mass value

    Returns
    -------
    str
        The matplotlib marker string for the given mass

    Raises
    ------
    ValueError
        If the mass is not in the standard set of masses
    """
    # Round to nearest integer to handle floating point precision
    mass_rounded = round(mass)

    # Find the closest standard mass
    closest_mass = min(STANDARD_MASSES, key=lambda x: abs(x - mass_rounded))

    # Check if it's close enough (within 1% tolerance)
    if abs(closest_mass - mass_rounded) / closest_mass > 0.01:
        raise ValueError(
            f"Mass {mass} is not close to any standard mass: {STANDARD_MASSES}"
        )

    return MASS_MARKERS[closest_mass]


def get_mass_style(mass):
    """
    Get both color and marker for a given molecular mass.

    Parameters
    ----------
    mass : float
        The molecular mass value

    Returns
    -------
    tuple
        (color, marker) tuple for the given mass
    """
    return get_mass_color(mass), get_mass_marker(mass)


def scatter_with_mass_style(ax, x, y, mass, label=None, **kwargs):
    """
    Create a scatter plot with standardized color and marker for the given mass.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes to plot on
    x, y : array-like
        The x and y coordinates
    mass : float
        The molecular mass to determine color and marker
    label : str, optional
        The label for the scatter plot
    **kwargs
        Additional keyword arguments passed to ax.scatter

    Returns
    -------
    matplotlib.collections.PathCollection
        The scatter plot object
    """
    color, marker = get_mass_style(mass)

    # Use the provided label or create a default one
    if label is None:
        label = f"${mass:.0f}$ g/mol"

    return ax.scatter(x, y, c=[color], marker=marker, label=label, **kwargs)


def get_all_mass_legend_entries():
    """
    Get legend entries for all standard masses.

    Returns
    -------
    list
        List of tuples (mass, color, marker) for all standard masses
    """
    return [(mass, MASS_COLORS[mass], MASS_MARKERS[mass]) for mass in STANDARD_MASSES]


def format_mass_label(mass):
    """
    Format a mass value for use in labels and legends.

    Parameters
    ----------
    mass : float
        The molecular mass value

    Returns
    -------
    str
        Formatted label string
    """
    return f"${mass:.0f}$ g/mol"


def setup_ab_figure(shared_legend=True):
    """
    Setup a figure with two subplots for 'a)' and 'b)' style labeling.

    Parameters
    ----------
    shared_legend : bool, optional
        Whether to use a shared legend, default True

    Returns
    -------
    tuple
        (fig, axes) where axes is dict of subplot axes with keys 'A' and 'B' (and 'L' for the legend)
    """

    if shared_legend:
        # Calculate figure size for two square subplots
        # Height includes legend space (ratio 1:10) plus square subplot height
        subplot_size = 5  # Size for each square subplot
        fig_width = subplot_size * 2.0 * 1.0175  # Two subplots plus spacing
        fig_height = (
            subplot_size * 11 / 10
        )  # Square height plus legend space

        fig, axes = plt.subplot_mosaic(
            "LL;AB",
            figsize=(fig_width, fig_height),
            layout="constrained",
            height_ratios=[1, 10],
        )
    else:
        subplot_size = 5
        fig_width = subplot_size * 2.2
        fig_height = subplot_size

        fig, axes = plt.subplots(
            nrows=1, ncols=2, figsize=(fig_width, fig_height), layout="constrained"
        )

    fig.get_layout_engine().set(w_pad=0, h_pad=0, hspace=0.0175, wspace=0.0175)

    return fig, axes


def setup_figure():
    """
    Setup a basic figure with default size and layout.

    Returns
    -------
    tuple
        (fig, ax) where ax is the main axes object
    """
    fig, ax = plt.subplots(figsize=(7, 7), layout="constrained")
    fig.get_layout_engine().set(w_pad=0, h_pad=0, hspace=0.0175, wspace=0.0175)
    return fig, ax


def add_subplot_labels(axes, labels=None, position="top-left", offset=(0.02, 0.98)):
    """
    Add labels (like 'a)', 'b)') to subplots.

    Parameters
    ----------
    axes : list or array
        List of matplotlib axes objects
    labels : list, optional
        List of labels to use. If None, uses 'a)', 'b)', etc.
    position : str, optional
        Position of labels: 'top-left', 'top-right', 'bottom-left', 'bottom-right'
    offset : tuple, optional
        (x, y) offset in axes coordinates for label position
    """
    if labels is None:
        if isinstance(axes, dict):
            labels = [c.lower() + ")" for c in list(axes.keys())]
        else:
            labels = [f"{chr(97 + i)})" for i in range(len(axes))]

    # Define position mappings
    position_map = {
        "top-left": (offset[0], offset[1]),
        "top-right": (1 - offset[0], offset[1]),
        "bottom-left": (offset[0], offset[1] - 1),
        "bottom-right": (1 - offset[0], offset[1] - 1),
    }

    x_pos, y_pos = position_map.get(position, position_map["top-left"])

    for ax, label in zip(axes.values() if isinstance(axes, dict) else axes, labels):
        if label == "l)":
            # Special case for legend label
            continue
        # Add text label to the specified position in axes coordinates
        ax.text(
            x_pos,
            y_pos,
            label,
            transform=ax.transAxes,
            fontweight="bold",
            va="top",
            ha=position.split("-")[1],  # 'left' or 'right'
        )


def setup_shared_legend(fig, axes, ncol=None):
    """
    Setup a shared legend for all subplots.

    Parameters
    ----------
    axes : dict
        Matplotlib axes dictionary with key 'L' for legend
    ncol : int, optional
        Number of columns in legend. If None, finds a good number based on labels.
    """
    assert isinstance(axes, dict)
    assert "L" in axes, "Axes dictionary must contain a key 'L' for legend"

    # Collect all handles and labels from all axes
    handles, labels = [], []
    for ax in axes.values():
        h, label_list = ax.get_legend_handles_labels()
        for handle, label in zip(h, label_list):
            if label not in labels:  # Avoid duplicates
                handles.append(handle)
                labels.append(label)

    # Remove individual legends from axes
    for ax in axes.values():
        legend = ax.get_legend()
        if legend:
            legend.remove()

    # Set ncol if not provided
    if ncol is None:
        ncol = len(labels)
        if ncol > 4:
            # find a good multiple of 2 or 3
            ncol = (
                4 if ncol % 4 == 0 else 3 if ncol % 3 == 0 or (ncol + 1) % 3 == 0 else 2
            )

    # Add shared legend to figure
    axes["L"].legend(
        handles,
        labels,
        loc="center",
        mode="expand",
        ncol=ncol,
        frameon=True,
        borderaxespad=0.0,
        columnspacing=1.0,  # Reduce space between columns (default is 2.0)
        handletextpad=0.4,
        # Reduce space between symbol and text (default is 0.8)
    )

    # Make legend axes frame invisible
    axes["L"].axis("off")
