#!/usr/bin/env python
import os
import random
import warnings
from math import sqrt

import matplotlib.pyplot as plt
import numpy as np

# This script is more complicated than necessary,
# as it is originating from the hexagonal one
#
# we first generate the positions of the crosslinkers
# (the ones that actually form the hexagonal lattice)

chain_length = 51 #5**2  # 5 # 5**2
nr_of_cells_in_dir_x = 5  # 4
nr_of_cells_in_dir_y = 5  # 2
plot_boundary = True  # False
include_ghost_bonds = False
mc = False
random.seed(8803) # 555
ideal_functionality = 4

# a: side length, b: longer tilted side length, c: shorter tilted side length
a = chain_length

box_len = (nr_of_cells_in_dir_x)*a
box_height = (nr_of_cells_in_dir_y)*a

# setup the crosslinker coordinates
crosslinker_coordinates = []
bonds = []
current_y = 0
for ycell in range(int(nr_of_cells_in_dir_y)):
    current_x = 0
    for xcell in range(int(nr_of_cells_in_dir_x)):
        crosslinker_coordinates.append([
            current_x, current_y
        ])
        # next it
        current_x += a
    current_y += a


def compute_distance_between_coords(coords1, coords2):
    is_pbc = False
    deltax = abs(coords2[0] - coords1[0])
    if (deltax > 0.5*box_len):
        deltax -= box_len
        is_pbc = True
    deltay = abs(coords2[1] - coords1[1])
    if (deltay > 0.5*box_height):
        deltay -= box_height
        is_pbc = True
    return np.sqrt(deltax**2 + deltay**2), is_pbc


# setup the bonds
# they could easily be setup in an algorithm of O(N),
# yet the following is O(N^2). As it only has to run once,
# and is fast enough for our purposes, this comment
# shall serve as a warning if you want to use the algorithm
# for other purposes
available_sites = np.zeros(len(crosslinker_coordinates)) + ideal_functionality
max_dist = chain_length if not mc else chain_length*sqrt(chain_length)
for i in range(len(crosslinker_coordinates)):
    if (available_sites[i] <= 0):
        continue
    potential_bonds = []
    indices_to_compare = list(range(0, len(crosslinker_coordinates)))
    random.shuffle(indices_to_compare)
    for j in indices_to_compare:
        if (available_sites[j] <= 0):
            continue
        if (available_sites[i] <= 0):
            break
        if (i == j and available_sites[i] < 2):
            continue
        if (not mc and i == j):
            continue
        delta, is_periodic = compute_distance_between_coords(
            crosslinker_coordinates[i], crosslinker_coordinates[j])
        if (delta < 1.005*max_dist):
            sorted_res = sorted([i, j])
            available_sites[i] -= 1
            available_sites[j] -= 1
            sorted_res.append(is_periodic)
            bonds.append(tuple(sorted_res))

if (np.sum(available_sites) != 0.0):
    warnings.warn("Remaining available sites: {}".format(
        np.sum(available_sites)))
# remove duplicate bonds (every one should be duplicate currently)
# bonds = list(set(bonds))
assert(len(bonds) == 2 * nr_of_cells_in_dir_x * nr_of_cells_in_dir_y)

fig, ax = plt.subplots()
# fig.set_size_inches(boxLen/chainLength*0.875*1.5, boxHeight /
#                     chainLength*0.875*1.5, forward=True)
fig.set_size_inches(7.0, 6.062177826491079, forward=True)

crosslinker_coordinates = np.array(crosslinker_coordinates)

# the offset of the box we draw; bonds cannot be split in same lengths
box_offset_y = -(0.0)*a
box_offset_x = -(0.25)*a
box_coords_y = [box_offset_y, box_height+box_offset_y]
box_coords_x = [box_offset_x, box_len+box_offset_x]


def compute_periodic_atom_neighbour(atom_from, atom_to):
    new_atom = atom_to.copy()
    delta_x = abs(atom_from[0] - atom_to[0])
    if (delta_x > 0.5*box_len):
        if (atom_to[0] > 0.5*box_len):
            new_atom[0] -= box_len
        else:
            new_atom[0] += box_len
    delta_y = abs(atom_from[1] - atom_to[1])
    if (delta_y > 0.5*box_height):
        if (atom_to[1] > 0.5*box_height):
            new_atom[1] -= box_height
        else:
            new_atom[1] += box_height
    return new_atom


# plot cross-linker atoms
ax.scatter(crosslinker_coordinates[:, 0],
           crosslinker_coordinates[:, 1], s=89, zorder=2)

# plot box
line_thickness = 3
ax.plot([box_coords_x[0], box_coords_x[0]], [box_coords_y[0], box_coords_y[1]],  '-',
        color='tab:orange', linewidth=line_thickness, zorder=5)
ax.plot([box_coords_x[0], box_coords_x[1]], [box_coords_y[1], box_coords_y[1]], '-',
        color='tab:orange', linewidth=line_thickness, zorder=5)
ax.plot([box_coords_x[1], box_coords_x[1]], [box_coords_y[0], box_coords_y[1]], '-',
        color='tab:orange', linewidth=line_thickness, zorder=5)
ax.plot([box_coords_x[0], box_coords_x[1]], [box_coords_y[0], box_coords_y[0]], '-',
        color='tab:orange', linewidth=line_thickness, zorder=5)

# count functionality for verification purposes
crosslinker_functionality = {}
# plot bonds (and count functionality)
for i in range(len(bonds)):
    bond = bonds[i]
    if (bond[0] not in crosslinker_functionality):
        crosslinker_functionality[bond[0]] = 0
    crosslinker_functionality[bond[0]] += 1

    if (bond[1] not in crosslinker_functionality):
        crosslinker_functionality[bond[1]] = 0
    crosslinker_functionality[bond[1]] += 1

    atom_from = crosslinker_coordinates[bond[0]]
    atom_to = crosslinker_coordinates[bond[1]]
    if (not bond[2]):
        ax.plot([atom_from[0], atom_to[0]], [atom_from[1], atom_to[1]], '-',
                color="black", zorder=1)
    else:
        # plot the "cut-off", periodic bonds
        def plot_line_half_dashed(afrom, ato, color="black"):
            x0 = afrom[0]
            x1 = ato[0]
            if (x1 == x0):
                x1 = x0 + 0.5*a
            y0 = afrom[1]
            y1 = ato[1]
            if (y1 == y0):
                y1 = x0 + 0.5*a
            slope = (y0-y1)/(x0-x1)
            intercept = y1 - slope*x1

            def y(x):
                return slope*x + intercept

            def x(y):
                if (slope == 0):
                    return x1
                return (y-intercept)/slope

            def find_nearest_point(origin, target, points_to_check):
                nearest_point = None
                nearest_point_d = None
                for point in points_to_check:
                    # TODO: limit: only within box spawned by x0, y0, x1, y1
                    if (point[0] > origin[0] and point[0] > target[0]):
                        continue
                    if (point[0] < origin[0] and point[0] < target[0]):
                        continue
                    if (point[1] > origin[1] and point[1] > target[1]):
                        continue
                    if (point[1] < origin[1] and point[1] < target[1]):
                        continue

                    current_point_d = (point[0]-origin[0])**2 + (point[1]-origin[1])**2
                    if (nearest_point_d is None or current_point_d < nearest_point_d):
                        nearest_point_d = current_point_d
                        nearest_point = point

                if (nearest_point is None):
                    nearest_point = target

                return nearest_point

            # find the actually relevant point
            point_of_intercept = find_nearest_point((x0, y0), (x1, y1), [
                (box_coords_x[0], y(box_coords_x[0])),
                (box_coords_x[1], y(box_coords_x[1])),
                (x(box_coords_y[0]), box_coords_y[0]),
                (x(box_coords_y[1]), box_coords_y[1])])

            # solid part
            ax.plot([x0, point_of_intercept[0]], [
                y0, point_of_intercept[1]], '-', color=color, zorder=1)
            # dashed part
            ax.plot([point_of_intercept[0], x1], [
                point_of_intercept[1], y1], ':', color=color, zorder=1)

        # need to plot twice the same bond, half
        atom_new_to = compute_periodic_atom_neighbour(atom_from, atom_to)
        atom_new_from = compute_periodic_atom_neighbour(atom_to, atom_from)

        plot_line_half_dashed(atom_from, atom_new_to)
        plot_line_half_dashed(atom_to, atom_new_from)

# for i in range(len(crosslinkerCoordinates)):
#     ax.annotate(str(i), crosslinkerCoordinates[i])

if (not plot_boundary):
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    ax.set_frame_on(False)

# ax.set(title="len: {}, nX: {}, nY: {}".format(
#     chainLength, nr_of_cells_in_dir_x, nr_of_cells_in_dir_y))
fig.savefig(os.path.join(os.path.dirname(__file__), "plots",
            "figure{}-2d-square-lattice-{}{}-{}x{}.png".format("-mc" if mc else "", "no-boundary-" if not plot_boundary else "", chain_length, nr_of_cells_in_dir_x, nr_of_cells_in_dir_y)), bbox_inches="tight", dpi=300)

# verify the functionality of each atom
for i in range(len(crosslinker_coordinates)):
    if (crosslinker_functionality[i] != ideal_functionality):
        warnings.warn("Found crosslinker {} with functionality {} != {}".format(
            i, crosslinker_functionality[i], ideal_functionality))

# add "in-between" atoms on every bond
actual_bonds = []
other_atom_coordinates = []
basis_id = len(crosslinker_coordinates)
for i in range(len(bonds)):
    bond = bonds[i]
    atom_from = crosslinker_coordinates[bond[0]]
    atom_to = crosslinker_coordinates[bond[1]]
    if (bond[2]):
        # periodic case: we cheat a bit. Let the atoms escape the box, they will be wrapped by LAMMPS
        atom_to = compute_periodic_atom_neighbour(atom_from, atom_to)
    # non-periodic case: "easy"
    delta_x = atom_to[0] - atom_from[0]
    delta_y = atom_to[1] - atom_from[1]
    current_x = atom_from[0]
    current_y = atom_from[1]
    for new_atom_i in range(chain_length-1):
        current_x += delta_x/chain_length
        current_y += delta_y/chain_length
        other_atom_coordinates.append([
            current_x, current_y
        ])
        this_atom_id = basis_id+len(other_atom_coordinates)-1
        if (new_atom_i == 0):
            # opening one
            actual_bonds.append([bond[0], this_atom_id])
        else:
            actual_bonds.append([this_atom_id-1, this_atom_id])
    # and the closing one
    actual_bonds.append([bond[1], basis_id+len(other_atom_coordinates)-1])

# finally, write the LAMMPS data file
print("LAMMPS 2D square lattice data file: nx = {}, ny = {}, chain len = {}".format(
    nr_of_cells_in_dir_x, nr_of_cells_in_dir_y, chain_length))
print("""
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
""".format(len(other_atom_coordinates)+len(crosslinker_coordinates), 
           len(actual_bonds)+(len(bonds) if include_ghost_bonds else 0), 2 if include_ghost_bonds else 1))
print("""
      {}       {} xlo xhi
      {}       {} ylo yhi
      {}       {} zlo zhi
""".format(0, box_len, 0, box_height, -0.5, 0.5))
print("""
Masses

  1  1.0
  2  1.0

Atoms
""")
molecule_idx = 1
atom_idx = 1
for i in range(len(crosslinker_coordinates)):
    print("         {}       {}       {}    {}    {}    {}   0   0   0".format(
        atom_idx, molecule_idx, 2, crosslinker_coordinates[i][0], crosslinker_coordinates[i][1], 0))
    atom_idx += 1

for i in range(len(other_atom_coordinates)):
    print("         {}       {}       {}    {}    {}    {}   0   0   0".format(
        atom_idx, molecule_idx, 1, other_atom_coordinates[i][0], other_atom_coordinates[i][1], 0))
    atom_idx += 1

print("""
Bonds
""")

bond_idx = 1
for i in range(len(actual_bonds)):
    print("        {}  1        {}        {}".format(
        i, actual_bonds[i][0]+1, actual_bonds[i][1]+1))
    bond_idx += 1

# also print the bonds of the crosslinkers as zero bonds
if (include_ghost_bonds):
    for i in range(len(bonds)):
        print("        {}  2        {}        {}".format(
            i, bonds[i][0]+1, bonds[i][1]+1))
        bond_idx += 1


if (include_ghost_bonds):
    print("""
    Bond Coeffs

    1 harmonic 4.1419333334 0
    """)
    print("2 zero")
