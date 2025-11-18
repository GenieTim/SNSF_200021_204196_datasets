#!/usr/bin/env python

import math
import os
import random
import warnings

import numpy as np
import pandas as pd
from numpy.random import default_rng
# https://polymercpp.readthedocs.io/en/latest/
from PolymerCpp.helpers import getCppSAWLC
from pylimer_tools.calc.structure_analysis import \
    compute_stoichiometric_imbalance
from pylimer_tools.io.read_lammps_output_file import read_data_file
from pylimer_tools_cpp import Box, DataFileWriter, Universe

# input
zero_coords_only = True

# main chains
main_beads_per_chain = 25  # 80
nr_of_chains = 500  # 100
n_crosslinks = 0 # int((2*100)/2)  # int(round(1.61*100/(2)))  # int(round(210/2*1.71))
bead_distance = 1.0
regular_atom_type = 1
crosslink_type = 2
use_type_per_molecule = False
# 1.23  # 1.9206  # None if you do not want a distribution
polydispersity_index = None
density = 0.85  # beads/sigma^3
versions = 1
dimensionality = 3

# monofunctional chains
beads_per_monofunctional_chain = 442
nr_of_monofunctional_chains = 0  # 48
monofunctional_bead_type = 5

# solvent chains
beads_per_solvent_chain = 19
nr_of_solvent_chains = 0
solvent_atom_type = 3

#
assert (dimensionality == 2 or dimensionality == 3)

for version in range(1, versions + 1):
    seed = 8804 + main_beads_per_chain + nr_of_chains + version
    # whether to number the molecules per chain or per position in chain
    mode_swappable = False
    file_to_write = os.path.join(os.path.dirname(__file__), "structure", "melt_{}_a_{}{}{}{}{}{}{}_v_{}.{}d_structure.out".format(
        nr_of_chains, main_beads_per_chain,
        "_{}_xlinks".format(n_crosslinks) if n_crosslinks > 0 else "",
        "_I_{}".format(
            polydispersity_index) if polydispersity_index is not None else "",
        "_solv_{}_a_{}".format(
            nr_of_solvent_chains, beads_per_solvent_chain) if nr_of_solvent_chains > 0 else "",
        "_mono_{}_a_{}".format(nr_of_monofunctional_chains,
                               beads_per_monofunctional_chain) if nr_of_monofunctional_chains > 0 else "",
        "_type_per_molecule" if use_type_per_molecule else "",
        "_zero_coords" if zero_coords_only else "",
        version, dimensionality)
    )

    # actual scripts
    random.seed(seed)
    if (nr_of_solvent_chains > 0 and use_type_per_molecule):
        warnings.warn(
            "CAUTION: it might be that your solvent is inconsistent with the normal beads")

    nr_of_beads_per_chain = np.repeat(main_beads_per_chain, nr_of_chains)
    if (polydispersity_index is not None):
        sigma = np.sqrt(np.log(polydispersity_index))
        mu = np.log(main_beads_per_chain) - (sigma**2) / 2
        rng = default_rng(seed)
        nr_of_beads_per_chain = rng.lognormal(
            mu, sigma, nr_of_chains)
        nr_of_beads_per_chain = np.array([int(round(i))
                                          for i in nr_of_beads_per_chain])
        print("Mean nr of beads per chain: {} instead of {}".format(
            np.mean(nr_of_beads_per_chain), main_beads_per_chain))

    total_nr_of_beads = np.sum(nr_of_beads_per_chain) + beads_per_solvent_chain * nr_of_solvent_chains + \
        n_crosslinks + nr_of_monofunctional_chains * beads_per_monofunctional_chain
    target_side_len = ((total_nr_of_beads) / density)**(1. / 3.)
    base_side_len = target_side_len * ((1 / bead_distance))
    universe = Universe(base_side_len, base_side_len, base_side_len)

    origin_displacement = nr_of_chains / main_beads_per_chain
    current_max_id = 1

    # assemble all to add at once
    atoms_to_add = {
        "ids": [], "types": [], "x": [], "y": [], "z": []
    }
    bonds_to_add = {
        "from": [], "to": []
    }
    angles_to_add = {
        "from": [], "to": [], "via": []
    }

    def add_chain_to_universe(num_beads_in_chain, bead_type, final_bead_type=None):
        # do actual random walking
        chain = getCppSAWLC(num_beads_in_chain - 1, 0.5 * bead_distance, 0.75)
        assert (len(chain) == num_beads_in_chain)
        chain_df = pd.DataFrame(chain, columns=list('xyz'))
        chain_df[list('xyz')] += (i - origin_displacement) * \
            np.array([random.random(), random.random(), random.random()])
        atom_ids = np.array(range(len(chain))) + current_max_id

        atoms_to_add['ids'].extend(atom_ids)
        atoms_to_add['types'].extend(np.repeat(bead_type, len(chain)))
        if (final_bead_type is not None):
            atoms_to_add["types"][-1] = final_bead_type
        atoms_to_add['x'].extend(chain_df['x'] if not zero_coords_only else [0. for _ in range(len(atom_ids))])
        atoms_to_add['y'].extend(chain_df['y'] if not zero_coords_only else [0. for _ in range(len(atom_ids))])
        atoms_to_add['z'].extend(chain_df['z'] if (dimensionality == 3 and not zero_coords_only) else 
                                 [0. for _ in range(len(atom_ids))])

        bonds_to_add['from'].extend(atom_ids[0:-1])
        bonds_to_add['to'].extend(atom_ids[1:])

        angles_to_add['from'].extend(atom_ids[0:-2])
        angles_to_add['to'].extend(atom_ids[1:-1])
        angles_to_add['via'].extend(atom_ids[2:])

    masses = {}
    print("Adding ordinary chains...")

    for i in range(nr_of_chains):
        atom_type = regular_atom_type if not use_type_per_molecule else i + regular_atom_type
        add_chain_to_universe(
            nr_of_beads_per_chain[i], atom_type)
        if (atom_type not in masses):
            masses[atom_type] = 1.
        current_max_id += nr_of_beads_per_chain[i]

    if (nr_of_solvent_chains > 0):
        print("Adding solvent chains...")
        masses[solvent_atom_type] = 1.

        for i in range(nr_of_solvent_chains):
            add_chain_to_universe(beads_per_solvent_chain, solvent_atom_type)
            current_max_id += beads_per_solvent_chain

    if (nr_of_monofunctional_chains > 0):
        print("Adding monofunctional chains")
        masses[monofunctional_bead_type] = 1.

        for i in range(nr_of_monofunctional_chains):
            add_chain_to_universe(beads_per_monofunctional_chain,
                                  regular_atom_type, monofunctional_bead_type)
            current_max_id += beads_per_monofunctional_chain

    if (n_crosslinks > 0):
        print("Adding cross-links")
        masses[crosslink_type] = 1.
        for i in range(n_crosslinks):
            atom_id = 1 + current_max_id

            atoms_to_add['ids'].append(atom_id)
            atoms_to_add['types'].append(crosslink_type)
            atoms_to_add['x'].append(random.random() * base_side_len)
            atoms_to_add['y'].append(random.random() * base_side_len)
            atoms_to_add['z'].append(random.random() * base_side_len)

            current_max_id += 1

    # adjust masses for LAMMPS
    for i in range(1, max([crosslink_type if n_crosslinks > 0 else 0,
                           monofunctional_bead_type if nr_of_monofunctional_chains > 0 else 0,
                           regular_atom_type, solvent_atom_type if nr_of_solvent_chains > 0 else 0])):
        if i not in masses:
            masses[i] = 1.

    print("Assembling universe...")
    universe.set_masses(masses)

    # we do not care about the image -> our coordinates escape
    zeros = np.repeat(0, len(atoms_to_add['ids']))
    universe.add_atoms(atoms_to_add['ids'], atoms_to_add['types'], atoms_to_add['x'],
                       atoms_to_add['y'], atoms_to_add['z'], zeros, zeros, zeros)
    universe.add_bonds(bonds_to_add['from'], bonds_to_add['to'])
    universe.add_angles(angles_to_add['from'],
                        angles_to_add['via'], angles_to_add['to'], np.repeat(1, len(angles_to_add['from'])))

    # achieve the target bond distance
    universe.set_box(Box(target_side_len, target_side_len,
                         target_side_len), rescale_atoms=True)

    # verify system
    print("Verifying system...")
    assert (universe.get_nr_of_angles() == np.sum(nr_of_beads_per_chain - 2) +
            (nr_of_solvent_chains * (beads_per_solvent_chain - 2)) +
            (nr_of_monofunctional_chains * (beads_per_monofunctional_chain - 2)))

    mean_bond_len = np.mean([np.mean(m.compute_bond_lengths())
                            for m in universe.get_molecules(2)])
    print("Mean bond length: {}".format(mean_bond_len))

    # output new system
    data_writer = DataFileWriter(universe)
    data_writer.config_include_angles(True)
    data_writer.config_include_dihedral_angles(False)
    data_writer.config_molecule_idx_for_swap(mode_swappable)
    data_writer.config_move_into_box(True)
    data_writer.config_attempt_image_reset(True)
    data_writer.write_to_file(file_to_write)
    print("Wrote file \"{}\"".format(file_to_write))

    print("Verifying file...")
    new_universe = read_data_file(file_to_write)
    assert (universe.get_nr_of_atoms() == new_universe.get_nr_of_atoms())
    assert (len(universe.get_atoms_by_type(2)) ==
            len(new_universe.get_atoms_by_type(2)))
    assert (1e-3 > abs(mean_bond_len -
                       np.mean([np.mean(m.compute_bond_lengths()) for m in new_universe.get_molecules(2)])))
    assert (math.isclose(new_universe.get_nr_of_atoms() /
                         new_universe.get_volume(), density))
    print("Density is {}".format(
        new_universe.get_nr_of_atoms() / new_universe.get_volume()))
    if (polydispersity_index is not None):
        print("Polydispersity is {}/{}".format(
            new_universe.compute_polydispersity_index(crosslink_type), polydispersity_index))
    if (n_crosslinks > 0):
        print("r = {}".format(compute_stoichiometric_imbalance(new_universe, crosslink_type, functionality_per_type={
            regular_atom_type: 2,
            crosslink_type: 4,
            monofunctional_bead_type: 1
        })))
    # TODO: some more
    print("Verification passed.")
