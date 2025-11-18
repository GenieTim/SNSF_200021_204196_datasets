#!/usr/bin/env python
from math import sqrt
import os
import warnings
import random

import matplotlib.pyplot as plt
import numpy as np

# we first generate the positions of the crosslinkers
# (the ones that actually form the hexagonal lattice)

chain_length = 26  # 5 # 5**2
nr_of_cells_in_dir_x = 5
nr_of_cells_in_dir_y = 5
plot_boundary = False
include_ghost_bonds = False
random.seed(5)  # 3

# a: side length, b: longer tilted side length, c: shorter tilted side length
a = chain_length
b = np.sqrt(3) / 2 * a
c = 0.5 * a

# make sure the periodicity can be achieved
if nr_of_cells_in_dir_x % 2 != 0:
    nr_of_cells_in_dir_x += 1
if nr_of_cells_in_dir_y % 2 != 0:
    nr_of_cells_in_dir_y += 1

box_len = (nr_of_cells_in_dir_x * 2) * b
box_height = (nr_of_cells_in_dir_y) * (a + c)

# setup the crosslinker coordinates
crosslinker_coordinates = []
bonds = []
current_y = a
for ycell in range(int(nr_of_cells_in_dir_y / 2)):
    current_x = 0.5 * b
    for xcell in range(int(nr_of_cells_in_dir_x)):
        crosslinker_coordinates.append([current_x, current_y])
        crosslinker_coordinates.append([current_x + b, current_y - c])
        crosslinker_coordinates.append([current_x, current_y + a])
        crosslinker_coordinates.append([current_x + b, current_y + a + c])
        current_x += 2 * b
    current_y += 2 * (c + a)


def compute_distance_between_coords(coords1, coords2):
    is_pbc = False
    deltax = abs(coords2[0] - coords1[0])
    if deltax > 0.5 * box_len:
        deltax -= box_len
        is_pbc = True
    deltay = abs(coords2[1] - coords1[1])
    if deltay > 0.5 * box_height:
        deltay -= box_height
        is_pbc = True
    return np.sqrt(deltax**2 + deltay**2), is_pbc


def compute_periodic_atom_neighbour(atom_from, atom_to):
    new_atom = atom_to.copy()
    deltax = abs(atom_from[0] - atom_to[0])
    if deltax > 0.5 * box_len:
        if atom_to[0] > 0.5 * box_len + 0.5 * b:
            new_atom[0] -= box_len
        else:
            new_atom[0] += box_len
    deltay = abs(atom_from[1] - atom_to[1])
    if deltay > 0.5 * box_height:
        if atom_to[1] > 0.5 * box_height + a:
            new_atom[1] -= box_height
        else:
            new_atom[1] += box_height
    return new_atom


# setup the bonds
# they could easily be setup in an algorithm of O(N),
# yet the following is O(N^2). As it only has to run once,
# and is fast enough for our purposes, this comment
# shall serve as a warning if you want to use the algorithm
# for other purposes
available_sites = np.zeros(len(crosslinker_coordinates)) + 3
for i in range(len(crosslinker_coordinates)):
    if available_sites[i] <= 0:
        continue
    potential_bonds = []
    indices_to_compare = list(range(0, len(crosslinker_coordinates)))
    random.shuffle(indices_to_compare)
    for j in indices_to_compare:
        if available_sites[i] <= 0:
            break
        if available_sites[j] <= 0:
            continue
        if i == j and available_sites[i] < 2:
            continue
        delta, is_periodic = compute_distance_between_coords(
            crosslinker_coordinates[i], crosslinker_coordinates[j]
        )
        if delta < 1.005 * chain_length * sqrt(chain_length):
            sorted_res = sorted([i, j])
            available_sites[i] -= 1
            available_sites[j] -= 1
            sorted_res.append(is_periodic)
            bonds.append(tuple(sorted_res))

# remove duplicate bonds (every one should be duplicate currently)
bonds = list(set(bonds))

fig, ax = plt.subplots()
fig.set_size_inches(
    box_len / chain_length * 0.875, box_height / chain_length * 0.875, forward=True
)

crosslinker_coordinates = np.array(crosslinker_coordinates)

# count functionality for verification purposes
crosslinker_functionality = {}
# plot bonds (and count functionality)
for i in range(len(bonds)):
    bond = bonds[i]
    if bond[0] not in crosslinker_functionality:
        crosslinker_functionality[bond[0]] = 0
    crosslinker_functionality[bond[0]] += 1

    if bond[1] not in crosslinker_functionality:
        crosslinker_functionality[bond[1]] = 0
    crosslinker_functionality[bond[1]] += 1

    atom_from = crosslinker_coordinates[bond[0]]
    atom_to = crosslinker_coordinates[bond[1]]
    if not bond[2]:
        ax.plot(
            [atom_from[0], atom_to[0]],
            [atom_from[1], atom_to[1]],
            "-",
            color="black",
            zorder=1,
        )
    else:
        # plot the "cut-off", periodic bonds
        def plot_line_half_dashed(afrom, ato, color="black"):
            x0 = afrom[0]
            x1 = ato[0]
            y0 = afrom[1]
            y1 = ato[1]
            ax.plot(
                [x0, x0 + (x1 - x0) / 2],
                [y0, y0 + (y1 - y0) / 2],
                "-",
                color=color,
                zorder=1,
            )
            # dashed part
            ax.plot(
                [x0 + (x1 - x0) / 2, x1],
                [y0 + (y1 - y0) / 2, y1],
                ":",
                color=color,
                zorder=1,
            )

        # need to plot twice the same bond, half
        atom_new_to = compute_periodic_atom_neighbour(atom_from, atom_to)
        atom_new_from = compute_periodic_atom_neighbour(atom_to, atom_from)

        plot_line_half_dashed(atom_from, atom_new_to)
        plot_line_half_dashed(atom_to, atom_new_from)

# plot atoms
ax.scatter(crosslinker_coordinates[:, 0], crosslinker_coordinates[:, 1], s=89, zorder=2)

# for i in range(len(crosslinkerCoordinates)):
#     ax.annotate(str(i), crosslinkerCoordinates[i])

# plot box
line_thickness = 3
ax.plot([0, 0], [0, box_height], "-", color="tab:orange", linewidth=line_thickness)
ax.plot(
    [0, box_len],
    [box_height, box_height],
    "-",
    color="tab:orange",
    linewidth=line_thickness,
)
ax.plot(
    [box_len, box_len],
    [0, box_height],
    "-",
    color="tab:orange",
    linewidth=line_thickness,
)
ax.plot([0, box_len], [0, 0], "-", color="tab:orange", linewidth=line_thickness)
if not plot_boundary:
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    ax.set_frame_on(False)

# plot unit cell
unit_cell_line_thickness = 2  # lineThickness/2
small_adjust = 0.0  # 0.25
ax.plot(
    [0 + small_adjust, 0 + small_adjust],
    [0 + small_adjust, 3 * a + small_adjust],
    "--",
    color="tab:green",
    linewidth=unit_cell_line_thickness,
)
ax.plot(
    [0 + small_adjust, 2 * b + small_adjust],
    [0 + small_adjust, 0 + small_adjust],
    "--",
    color="tab:green",
    linewidth=unit_cell_line_thickness,
)
ax.plot(
    [2 * b + small_adjust, 2 * b + small_adjust],
    [0 + small_adjust, 3 * a + small_adjust],
    "--",
    color="tab:green",
    linewidth=unit_cell_line_thickness,
)
ax.plot(
    [0 + small_adjust, 2 * b + small_adjust],
    [3 * a + small_adjust, 3 * a + small_adjust],
    "--",
    color="tab:green",
    linewidth=unit_cell_line_thickness,
)

# ax.set(title="len: {}, nX: {}, nY: {}".format(
#     chainLength, nrOfCellsInDirX, nrOfCellsInDirY))
fig.savefig(
    os.path.join(
        os.path.dirname(__file__),
        "plots",
        "figure-2d-hexagonal-mc-lattice-{}{}-{}x{}.png".format(
            "no-boundary-" if not plot_boundary else "",
            chain_length,
            nr_of_cells_in_dir_x,
            nr_of_cells_in_dir_y,
        ),
    )
)

# verify the functionality of each atom
for i in range(len(crosslinker_coordinates)):
    if crosslinker_functionality[i] != 3:
        warnings.warn(
            "Found crosslinker {} with functionality {} != 3".format(
                i, crosslinker_functionality[i]
            )
        )

# add "in-between" atoms on every bond
actual_bonds = []
other_atom_coordinates = []
basis_id = len(crosslinker_coordinates)
for i in range(len(bonds)):
    bond = bonds[i]
    atom_from = crosslinker_coordinates[bond[0]]
    atom_to = crosslinker_coordinates[bond[1]]
    if bond[2]:
        # periodic case: we cheat a bit. Let the atoms escape the box, they will be wrapped by LAMMPS
        atom_to = compute_periodic_atom_neighbour(atom_from, atom_to)
    # non-periodic case: "easy"
    delta_x = atom_to[0] - atom_from[0]
    delta_y = atom_to[1] - atom_from[1]
    current_x = atom_from[0]
    current_y = atom_from[1]
    for new_atom_i in range(chain_length - 1):
        current_x += delta_x / chain_length
        current_y += delta_y / chain_length
        other_atom_coordinates.append([current_x, current_y])
        this_atom_id = basis_id + len(other_atom_coordinates) - 1
        if new_atom_i == 0:
            # opening one
            actual_bonds.append([bond[0], this_atom_id])
        else:
            actual_bonds.append([this_atom_id - 1, this_atom_id])
    # and the closing one
    actual_bonds.append([bond[1], basis_id + len(other_atom_coordinates) - 1])

# finally, write the LAMMPS data file
print(
    "LAMMPS hexagonal lattice data file: nx = {}, ny = {}, chain len = {}".format(
        nr_of_cells_in_dir_x, nr_of_cells_in_dir_y, chain_length
    )
)
print(
    """
          {} atoms
          {} bonds
          0  angles
          0  dihedrals
          0  impropers

          2  atom types
          {}  bond types
          0  angle types
          0  dihedral types
          0  improper types
""".format(
        len(other_atom_coordinates) + len(crosslinker_coordinates),
        len(actual_bonds) + (len(bonds) if include_ghost_bonds else 0),
        2 if include_ghost_bonds else 1,
    )
)
print(
    """
      {}       {} xlo xhi
      {}       {} ylo yhi
      {}       {} zlo zhi
""".format(
        0, box_len, 0, box_height, -0.5, 0.5
    )
)
print(
    """
Masses

  1  1.0
  2  1.0

Atoms
"""
)
molecule_idx = 1
atom_idx = 1
for i in range(len(crosslinker_coordinates)):
    print(
        "         {}       {}       {}    {}    {}    {}   0   0   0".format(
            atom_idx,
            molecule_idx,
            2,
            crosslinker_coordinates[i][0],
            crosslinker_coordinates[i][1],
            0,
        )
    )
    atom_idx += 1

for i in range(len(other_atom_coordinates)):
    print(
        "         {}       {}       {}    {}    {}    {}   0   0   0".format(
            atom_idx,
            molecule_idx,
            1,
            other_atom_coordinates[i][0],
            other_atom_coordinates[i][1],
            0,
        )
    )
    atom_idx += 1

print(
    """
Bonds
"""
)

bond_idx = 1
for i in range(len(actual_bonds)):
    print(
        "        {}  1        {}        {}".format(
            i, actual_bonds[i][0] + 1, actual_bonds[i][1] + 1
        )
    )
    bond_idx += 1

# also print the bonds of the crosslinkers as zero bonds
if include_ghost_bonds:
    for i in range(len(bonds)):
        print(
            "        {}  2        {}        {}".format(
                i, bonds[i][0] + 1, bonds[i][1] + 1
            )
        )
        bond_idx += 1


print(
    """
Bond Coeffs

1 harmonic 4.1419333334 0
"""
)
if include_ghost_bonds:
    print("2 zero")
