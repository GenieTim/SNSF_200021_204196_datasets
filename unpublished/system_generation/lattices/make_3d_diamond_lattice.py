#!/usr/bin/env python

import math
import os
import random
import sys
import warnings

import numpy as np
import pandas as pd
from numpy.random import default_rng

# https://polymercpp.readthedocs.io/en/latest/
from PolymerCpp.helpers import getCppSAWLC
from pylimer_tools.io.read_lammps_output_file import read_data_file
from pylimer_tools_cpp import (
    Box,
    DataFileWriter,
    MoleculeType,
    NeighbourList,
    Universe,
    do_random_walk_chain_from_to,
)

# input
target_density = 0.85

nr_of_cells_in_x = 5
nr_of_cells_in_y = nr_of_cells_in_x
nr_of_cells_in_z = nr_of_cells_in_x

nr_of_beads_per_chain_excluding_crosslinks = 25

crosslink_type = 2
regular_bead_type = 1

versions = 1

target_bond_len = 1.0

# Diamond cubic lattice
basic_diamond_cell_box = [4, 4, 4]
basic_diamond_cell_coords = [
    [0, 0, 0],
    [0, 2, 2],
    [2, 0, 2],
    [2, 2, 0],
    [3, 3, 3],
    [3, 1, 1],
    [1, 3, 1],
    [1, 1, 3],
]

nr_of_cells = nr_of_cells_in_x * nr_of_cells_in_y * nr_of_cells_in_z
nr_of_chains = nr_of_cells * 16
total_nr_of_atoms = (
    nr_of_chains * nr_of_beads_per_chain_excluding_crosslinks
    + nr_of_cells * len(basic_diamond_cell_coords)
)
target_side_len = (total_nr_of_atoms / target_density) ** (1.0 / 3.0)
scaling_factor = (
    np.mean(basic_diamond_cell_box)
    * np.mean([nr_of_cells_in_x, nr_of_cells_in_y, nr_of_cells_in_z])
    / target_side_len
)

for version in range(versions):
    atoms_to_add = {"ids": [], "types": [], "x": [], "y": [], "z": []}
    current_atom_id = 1
    for x in range(nr_of_cells_in_x):
        for y in range(nr_of_cells_in_y):
            for z in range(nr_of_cells_in_z):
                for i in range(len(basic_diamond_cell_coords)):
                    atoms_to_add["ids"].append(current_atom_id)
                    atoms_to_add["types"].append(crosslink_type)
                    atoms_to_add["x"].append(
                        basic_diamond_cell_coords[i][0] + x * basic_diamond_cell_box[0]
                    )
                    atoms_to_add["y"].append(
                        basic_diamond_cell_coords[i][1] + y * basic_diamond_cell_box[1]
                    )
                    atoms_to_add["z"].append(
                        basic_diamond_cell_coords[i][2] + z * basic_diamond_cell_box[2]
                    )
                    current_atom_id += 1

    xlink_only_universe = Universe(
        nr_of_cells_in_x * basic_diamond_cell_box[0],
        nr_of_cells_in_y * basic_diamond_cell_box[1],
        nr_of_cells_in_z * basic_diamond_cell_box[2],
    )
    zeros = np.repeat(0, len(atoms_to_add["ids"]))
    xlink_only_universe.add_atoms(
        atoms_to_add["ids"],
        atoms_to_add["types"],
        atoms_to_add["x"],
        atoms_to_add["y"],
        atoms_to_add["z"],
        zeros,
        zeros,
        zeros,
    )

    # we are lazy; while we analytically know where to set the chains,
    # we choose to use the neighbour list – faster to implement, less (visible) code
    # detect the bonds to make
    xlink_atoms = sorted(xlink_only_universe.get_atoms(), key=lambda a: a.get_id())
    neighbourlist = NeighbourList(
        xlink_atoms, xlink_only_universe.get_box(), math.sqrt(3.1)
    )

    chains_to_build = set()
    for xlink in xlink_atoms:
        neighbours = neighbourlist.get_atoms_close_to(xlink, math.sqrt(3.1), 0.0, False)
        assert len(neighbours) < 5
        for neigh in neighbours:
            if neigh.get_id() == xlink.get_id():
                continue
            chains_to_build.add(
                "{}-{}".format(
                    min(neigh.get_id(), xlink.get_id()),
                    max(neigh.get_id(), xlink.get_id()),
                )
            )

    assert len(chains_to_build) == 16 * nr_of_cells
    assert len(chains_to_build) == nr_of_chains

    # ok, now that we know which xlinks to connect, we can start interpolating between them
    bonds_to_add = {"from": [], "to": [], "type": []}

    # for now, we accept some random noise as sufficient to combat "ideal" equilibrium conditions
    for chain in chains_to_build:
        xlink_from = int(chain.split("-")[0]) - 1
        xlink_to = int(chain.split("-")[1]) - 1
        assert xlink_from != xlink_to

        vector_between_xlinks = xlink_atoms[xlink_from].compute_vector_to(
            xlink_atoms[xlink_to], xlink_only_universe.get_box()
        )
        random_walk_steps = do_random_walk_chain_from_to(
            xlink_only_universe.get_box(),
            [
                atoms_to_add["x"][xlink_from],
                atoms_to_add["y"][xlink_from],
                atoms_to_add["z"][xlink_from],
            ],
            [
                atoms_to_add["x"][xlink_to],
                atoms_to_add["y"][xlink_to],
                atoms_to_add["z"][xlink_to],
            ],
            nr_of_beads_per_chain_excluding_crosslinks,
            target_bond_len * scaling_factor,
        )
        assert len(random_walk_steps["x"]) == nr_of_beads_per_chain_excluding_crosslinks
        for i in range(nr_of_beads_per_chain_excluding_crosslinks):
            atoms_to_add["ids"].append(current_atom_id)
            atoms_to_add["types"].append(regular_bead_type)
            atoms_to_add["x"].append(random_walk_steps["x"][i])
            atoms_to_add["y"].append(random_walk_steps["y"][i])
            atoms_to_add["z"].append(random_walk_steps["z"][i])
            bonds_to_add["to"].append(current_atom_id)
            bonds_to_add["type"].append(1)
            if i == 0:
                bonds_to_add["from"].append(xlink_from + 1)
            else:
                bonds_to_add["from"].append(current_atom_id - 1)
            current_atom_id += 1
        bonds_to_add["to"].append(current_atom_id - 1)
        bonds_to_add["from"].append(xlink_to + 1)
        bonds_to_add["type"].append(1)

    universe = Universe(
        nr_of_cells_in_x * basic_diamond_cell_box[0],
        nr_of_cells_in_y * basic_diamond_cell_box[1],
        nr_of_cells_in_z * basic_diamond_cell_box[2],
    )
    universe.set_masses({regular_bead_type: 1, crosslink_type: 1})

    # we do not care about the image -> our coordinates escape
    assert len(atoms_to_add["ids"]) == total_nr_of_atoms
    zeros = np.repeat(0, len(atoms_to_add["ids"]))
    universe.add_atoms(
        atoms_to_add["ids"],
        atoms_to_add["types"],
        atoms_to_add["x"],
        atoms_to_add["y"],
        atoms_to_add["z"],
        zeros,
        zeros,
        zeros,
    )
    universe.add_bonds(bonds_to_add["from"], bonds_to_add["to"])

    angles = universe.detect_angles()
    universe.add_angles(
        angles["angle_from"],
        angles["angle_to"],
        angles["angle_via"],
        [1 for i in range(len(angles["angle_from"]))],
    )

    universe.set_box(
        Box(target_side_len, target_side_len, target_side_len), rescale_atoms=True
    )

    # verify system
    assert universe.get_nr_of_atoms() == total_nr_of_atoms
    mean_bond_len = np.mean(
        [
            np.mean(m.compute_bond_lengths())
            for m in universe.get_molecules(crosslink_type)
        ]
    )
    print("Mean bond length: {}".format(mean_bond_len))
    functionalities = universe.determine_effective_functionality_per_type()
    assert functionalities[regular_bead_type] == 2
    assert functionalities[crosslink_type] == 4

    molecules = universe.get_chains_with_crosslinker(crosslink_type)
    for m in molecules:
        assert m.get_strand_type() == MoleculeType.NETWORK_STRAND

    # output new system
    file_to_write = os.path.join(
        os.path.dirname(__file__),
        "structure",
        "3d-diamond-lattice_{}x{}x{}_a_{}_d_{}_v_{}.V-fixed.structure.out".format(
            nr_of_cells_in_x,
            nr_of_cells_in_y,
            nr_of_cells_in_z,
            nr_of_beads_per_chain_excluding_crosslinks,
            target_density,
            version,
        ),
    )
    data_writer = DataFileWriter(universe)
    data_writer.config_include_angles(True)
    data_writer.config_move_into_box(True)
    data_writer.config_attempt_image_reset(True)
    data_writer.write_to_file(file_to_write)
    print('Wrote file "{}"'.format(file_to_write))

    print("Verifying file...")
    new_universe = read_data_file(file_to_write)
    assert universe.get_nr_of_atoms() == new_universe.get_nr_of_atoms()
    assert len(universe.get_atoms_by_type(2)) == len(new_universe.get_atoms_by_type(2))
    assert len(new_universe.get_atoms_by_type(2)) == (
        nr_of_cells * len(basic_diamond_cell_coords)
    )
    assert [
        len(m) == nr_of_beads_per_chain_excluding_crosslinks
        for m in new_universe.get_molecules(2)
    ]
    assert len(new_universe.get_molecules(2)) == nr_of_chains
    assert 1e-3 > abs(
        mean_bond_len
        - np.mean(
            [np.mean(m.compute_bond_lengths()) for m in new_universe.get_molecules(2)]
        )
    )
    print(
        "Density is {}".format(
            new_universe.get_nr_of_atoms() / new_universe.get_volume()
        )
    )
    # TODO: some more
    print("Verification passed.")
