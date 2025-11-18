#!/usr/bin/env python

"""
"Butterfly" (closed starlike) structures: 
Structures of multiple connected link tuples, joined on one central link.
"""

import os
import random

import numpy as np
from pylimer_tools_cpp import DataFileWriter, MoleculeType, Universe

# input

# main chains
main_beads_per_chain = 25  # 80
nr_of_chains = 500  # 100

regular_atom_type = 1
crosslink_type = 2
nr_of_wings = 3

density = 0.85  # beads/sigma^3
versions = 1
dimensionality = 3

#
assert dimensionality == 2 or dimensionality == 3

for version in range(1, versions + 1):
    seed = 8804 + main_beads_per_chain + nr_of_chains + version
    # whether to number the molecules per chain or per position in chain
    mode_swappable = False
    file_to_write = os.path.join(
        os.path.dirname(__file__),
        "structure",
        "closed_starlike_{}_a_{}_f_{}_v_{}.{}d_structure.out".format(
            nr_of_chains, main_beads_per_chain, nr_of_wings * 2, version, dimensionality
        ),
    )

    # actual scripts
    random.seed(seed)

    # assemble all to add at once
    atoms_to_add = {"ids": [], "types": [], "x": [], "y": [], "z": []}
    bonds_to_add = {"from": [], "to": []}

    atom_id = 1
    for i in range(nr_of_chains):
        # start with cross-link
        atoms_to_add["ids"].append(atom_id)
        atoms_to_add["types"].append(crosslink_type)
        atoms_to_add["x"].append(0.25)
        atoms_to_add["y"].append(0.25)
        atoms_to_add["z"].append(0.25)

        last_crosslink_id = atom_id
        atom_id += 1
        for f in range(nr_of_wings):
            # add the spikes of the star
            for j in range(main_beads_per_chain * 2 + 1):
                atoms_to_add["ids"].append(atom_id)
                atoms_to_add["types"].append(
                    regular_atom_type if j != (main_beads_per_chain) else crosslink_type
                )
                atoms_to_add["x"].append(0.25)
                atoms_to_add["y"].append(0.25)
                atoms_to_add["z"].append(0.25)

                bonds_to_add["from"].append(
                    atom_id - 1 if j != 0 else last_crosslink_id
                )
                bonds_to_add["to"].append(atom_id)

                atom_id += 1
            # close the wing
            bonds_to_add["from"].append(atom_id - 1)
            bonds_to_add["to"].append(last_crosslink_id)
    # adjust masses for LAMMPS
    total_nr_of_beads = main_beads_per_chain * nr_of_chains + nr_of_chains
    target_side_len = ((total_nr_of_beads) / density) ** (1.0 / dimensionality)
    universe = Universe(
        target_side_len,
        target_side_len,
        target_side_len if dimensionality == 3 else 1.0,
    )

    masses = {}
    for i in range(1, max([crosslink_type, regular_atom_type])):
        if i not in masses:
            masses[i] = 1.0

    print("Assembling universe...")
    universe.set_masses(masses)

    # we do not care about the image -> our coordinates escape
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
    detected_angles = universe.detect_angles()
    universe.add_angles(
        detected_angles["angle_from"],
        detected_angles["angle_via"],
        detected_angles["angle_to"],
        [1 for _ in range(len(detected_angles["angle_from"]))],
    )

    mean_bond_len = np.mean(
        [np.mean(m.compute_bond_lengths()) for m in universe.get_molecules(2)]
    )
    print("Mean bond length: {}".format(mean_bond_len))

    functionalities = universe.determine_effective_functionality_per_type()
    assert functionalities[regular_atom_type] == 2
    # we have a central cross-link with f = 2 * nr_of_wings,
    # and nr_of_wings non-central cross-links with f = 2
    # for a total of nr_of_wings + 1 cross-links
    assert functionalities[crosslink_type] == (nr_of_wings * 4) / (nr_of_wings + 1)

    molecules = universe.get_chains_with_crosslinker(crosslink_type)
    for m in molecules:
        assert m.get_strand_type() == MoleculeType.NETWORK_STRAND
        assert m.get_nr_of_atoms() == main_beads_per_chain + 2

    # output new system
    data_writer = DataFileWriter(universe)
    data_writer.config_include_angles(True)
    data_writer.config_include_dihedral_angles(False)
    data_writer.config_molecule_idx_for_swap(mode_swappable)
    data_writer.config_move_into_box(True)
    data_writer.config_attempt_image_reset(True)
    data_writer.write_to_file(file_to_write)
    print('Wrote file "{}"'.format(file_to_write))
