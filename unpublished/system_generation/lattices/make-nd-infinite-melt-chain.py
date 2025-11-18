#!/usr/bin/env python

import os

import numpy as np

from pylimer_tools_cpp import (
    Box,
    DataFileWriter,
    MEHPForceRelaxation,
    MoleculeType,
    Universe,
)

# input
target_density = 0.85

nr_of_subsequent_chains_in_x = 50
nr_of_subsequent_chains_in_y = nr_of_subsequent_chains_in_x
nr_of_subsequent_chains_in_z = nr_of_subsequent_chains_in_x  # set to 0 for 2D

nr_of_parallel_chains = 5

nr_of_beads_per_chain_excluding_crosslinks = 25

crosslink_type = 2
regular_bead_type = 1

versions = 1

target_bond_len = 1.0

dimensionality = 3 if (nr_of_subsequent_chains_in_z > 0) else 2
assert nr_of_subsequent_chains_in_y == nr_of_subsequent_chains_in_x
assert (
    nr_of_subsequent_chains_in_z == 0
    or nr_of_subsequent_chains_in_z == nr_of_subsequent_chains_in_y
)

nr_of_chains = nr_of_parallel_chains * (
    nr_of_subsequent_chains_in_x
    + nr_of_subsequent_chains_in_y
    + nr_of_subsequent_chains_in_z
)
nr_of_chains = nr_of_chains * 1

total_nr_of_atoms = (
    nr_of_chains * nr_of_beads_per_chain_excluding_crosslinks + nr_of_chains
)

target_side_len = (total_nr_of_atoms / target_density) ** (1.0 / dimensionality)
base_side_len = nr_of_subsequent_chains_in_x * (
    nr_of_beads_per_chain_excluding_crosslinks + 1
)

for version in range(versions):
    atoms_to_add = {"ids": [], "types": [], "x": [], "y": [], "z": []}
    bonds_to_add = {"from": [], "to": []}
    current_atom_id = 1
    dirs = "xyz"
    n_chains_per_dir = [
        nr_of_subsequent_chains_in_x,
        nr_of_subsequent_chains_in_y,
        nr_of_subsequent_chains_in_z,
    ]
    for i in range(len(dirs)):
        # for each direction, make each chain
        for chain_in_dir_i in range(nr_of_parallel_chains):
            infinichain_start_atom_id = current_atom_id
            for subsequent_chain_i in range(n_chains_per_dir[i]):
                for atom_in_chain_i in range(
                    nr_of_beads_per_chain_excluding_crosslinks
                ):
                    atoms_to_add["ids"].append(current_atom_id)
                    atoms_to_add["types"].append(regular_bead_type)
                    for dir_i, dir in enumerate(dirs):
                        atoms_to_add[dir].append(
                            chain_in_dir_i * base_side_len / nr_of_parallel_chains
                            if dir_i != i
                            else atom_in_chain_i
                            + (
                                subsequent_chain_i
                                * (base_side_len)
                                / nr_of_subsequent_chains_in_x
                            )
                        )
                    bonds_to_add["from"].append(current_atom_id)
                    bonds_to_add["to"].append(current_atom_id + 1)
                    current_atom_id += 1
                atoms_to_add["ids"].append(current_atom_id)
                atoms_to_add["types"].append(crosslink_type)
                for dir_i, dir in enumerate(dirs):
                    atoms_to_add[dir].append(
                        chain_in_dir_i * base_side_len / nr_of_parallel_chains
                        if dir_i != i
                        else atom_in_chain_i
                        + (
                            subsequent_chain_i
                            * (base_side_len / nr_of_subsequent_chains_in_x)
                        )
                    )
                bonds_to_add["from"].append(current_atom_id)
                bonds_to_add["to"].append(
                    current_atom_id + 1
                    if subsequent_chain_i != n_chains_per_dir[i] - 1
                    else infinichain_start_atom_id
                )
                current_atom_id += 1
    # assembled

    universe = Universe(
        base_side_len, base_side_len, base_side_len if dimensionality == 3 else 1.0
    )
    universe.set_masses({regular_bead_type: 1, crosslink_type: 1})
    zeros = np.repeat(0, len(atoms_to_add["ids"]))
    universe.add_atoms(
        atoms_to_add["ids"],
        atoms_to_add["types"],
        atoms_to_add["x"],
        atoms_to_add["y"],
        atoms_to_add["z"] if dimensionality == 3 else zeros,
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
        Box(
            target_side_len,
            target_side_len,
            target_side_len if (dimensionality == 3) else 1.0,
        ),
        rescale_atoms=True,
    )

    functionalities = universe.determine_effective_functionality_per_type()
    assert functionalities[regular_bead_type] == 2
    assert functionalities[crosslink_type] == 2

    molecules = universe.get_chains_with_crosslinker(crosslink_type)
    for m in molecules:
        assert m.get_strand_type() == MoleculeType.NETWORK_STRAND
        assert m.get_nr_of_atoms() == nr_of_beads_per_chain_excluding_crosslinks + 2

    # output new system
    file_to_write = os.path.join(
        os.path.dirname(__file__),
        "structure",
        "infinite-melt-chain_{}x{}x{}_a_{}_v_{}.{}d_structure.out".format(
            nr_of_subsequent_chains_in_x,
            nr_of_subsequent_chains_in_y,
            nr_of_subsequent_chains_in_z,
            nr_of_beads_per_chain_excluding_crosslinks,
            version,
            dimensionality,
        ),
    )
    data_writer = DataFileWriter(universe)
    data_writer.config_include_angles(True)
    data_writer.config_move_into_box(True)
    data_writer.config_attempt_image_reset(True)
    data_writer.write_to_file(file_to_write)
    print('Wrote file "{}"'.format(file_to_write))

    force_relaxation = MEHPForceRelaxation(
        universe, remove_2functional_crosslinkers=False, remove_dangling_chains=False
    )
    force_relaxation.assume_box_large_enough(False)
    while force_relaxation.requires_another_run():
        force_relaxation.run_force_relaxation()
    # assert that we indeed have an infinite structure
    assert (
        force_relaxation.get_nr_of_active_springs()
        == force_relaxation.get_nr_of_springs()
    )
