#!/usr/bin/env python

import os

import numpy as np
# https://polymercpp.readthedocs.io/en/latest/
from pylimer_tools.io.read_lammps_output_file import read_data_file
from pylimer_tools_cpp import Box, DataFileWriter, Universe

# input
zero_coords_only = True

# main chains
main_beads_per_chain = 250  # 80
nr_of_tuples = 50  # 100
regular_atom_type = 1
crosslink_type = 2
bead_distance = 1
crosslink_functionality = 2
density = 0.5  # beads/sigma^3
versions = 1
dimensionality = 3

assert zero_coords_only
for version in range(1, versions + 1):

    file_to_write = os.path.join(
        os.path.dirname(__file__),
        "structure",
        "link_tuple_{}_a_{}_f_{}{}_v_{}.{}d_structure.out".format(
            nr_of_tuples,
            main_beads_per_chain,
            crosslink_functionality,
            "_zero_coords" if zero_coords_only else "",
            version,
            dimensionality,
        ),
    )

    # assemble all to add at once
    atoms_to_add = {"ids": [], "types": []}
    bonds_to_add = {"from": [], "to": []}

    current_atom_id = 1
    for tuple_idx in range(nr_of_tuples):
        atoms_to_add["ids"].append(current_atom_id)
        atoms_to_add["types"].append(crosslink_type)
        current_xlink1_id = current_atom_id
        current_atom_id += 1
        atoms_to_add["ids"].append(current_atom_id)
        atoms_to_add["types"].append(crosslink_type)
        current_xlink2_id = current_atom_id
        current_atom_id += 1
        for strand_idx in range(crosslink_functionality):
            # add a chain
            for atom_idx in range(main_beads_per_chain):
                atoms_to_add["ids"].append(current_atom_id)
                atoms_to_add["types"].append(regular_atom_type)

                if atom_idx == 0:
                    bonds_to_add["from"].append(current_xlink1_id)
                    bonds_to_add["to"].append(current_atom_id)
                else:
                    bonds_to_add["from"].append(current_atom_id - 1)
                    bonds_to_add["to"].append(current_atom_id)
                current_atom_id += 1
            bonds_to_add["from"].append(current_atom_id-1)
            bonds_to_add["to"].append(current_xlink2_id)
    total_nr_of_beads = current_atom_id - 1
    assert total_nr_of_beads == len(atoms_to_add["ids"])
    target_side_len = ((total_nr_of_beads) / density) ** (1.0 / dimensionality)
    universe = Universe(target_side_len, target_side_len, target_side_len)
    universe.set_box(
        Box(
            0,
            target_side_len,
            0,
            target_side_len,
            0 if dimensionality == 3 else -0.5,
            target_side_len if dimensionality == 3 else 0.5,
        )
    )

    zeros_i = np.repeat(0, len(atoms_to_add["ids"]))
    zeros_f = np.repeat(0.0, len(atoms_to_add["ids"]))
    # nor do we care about coordinates
    universe.add_atoms(
        atoms_to_add["ids"],
        atoms_to_add["types"],
        zeros_f,
        zeros_f,
        zeros_f,
        zeros_i,
        zeros_i,
        zeros_i,
    )
    universe.add_bonds(bonds_to_add["from"], bonds_to_add["to"])

    angles = universe.detect_angles()
    universe.add_angles(
        angles["angle_from"],
        angles["angle_to"],
        angles["angle_via"],
        [1 for i in range(len(angles["angle_from"]))],
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
    assert len(universe.get_atoms_by_type(crosslink_type)) == nr_of_tuples * 2
    assert (
        new_universe.determine_effective_functionality_per_type()[crosslink_type]
        == crosslink_functionality
    )
    assert (
        new_universe.determine_effective_functionality_per_type()[regular_atom_type]
        == 2
    )
    print("Verification passed.")
