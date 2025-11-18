#!/usr/bin/env python
import glob
import math
import os
import random

from pylimer_doctorate_utils.progressbar import progressbar
from pylimer_tools.io.read_lammps_output_file import read_data_file
from pylimer_tools_cpp import Atom, DataFileWriter

base_dir = os.path.dirname(__file__)

inputs = glob.glob(os.path.join(base_dir, "structure/*structure.out"))
random.shuffle(inputs)
crosslinker_type = 2
single_fixed_atom_type = 3
central_atom_type = 5
end_atom_type = 6

for input in progressbar(inputs):
    if "-cem" in input:
        continue
    universe = read_data_file(input)
    if len(universe.get_atoms_by_type(single_fixed_atom_type)):
        continue

    def update_atom_type(atom: Atom, new_type):
        global universe
        universe.replace_atom(
            atom.get_id(),
            Atom(
                atom.get_id(),
                new_type,
                atom.get_x(),
                atom.get_y(),
                atom.get_z(),
                atom.get_nx(),
                atom.get_ny(),
                atom.get_nz(),
            ),
        )
        assert universe.get_atom(atom.get_id()).get_type() == new_type
    n_crosslinks = len(universe.get_atoms_by_type(crosslinker_type))
    chains = universe.get_molecules(crosslinker_type)
    central_atoms_found = 0
    # mark "central" atoms and "end" atoms
    for chain in chains:
        atoms_lined_up = chain.get_atoms_lined_up()
        nr_of_atoms = len(atoms_lined_up)
        atoms_to_update = []
        if nr_of_atoms % 2 == 0:
            # need two center atoms
            atoms_to_update = [
                atoms_lined_up[int(nr_of_atoms / 2) - 1],
                atoms_lined_up[int(nr_of_atoms / 2)],
            ]
        else:
            atoms_to_update = [atoms_lined_up[math.floor(nr_of_atoms / 2)]]
        central_atoms_found += len(atoms_to_update)

        for atom in atoms_to_update:
            update_atom_type(atom, central_atom_type)
            assert universe.get_atom(atom.get_id()).get_type() == central_atom_type
        assert len(universe.get_atoms_by_type(central_atom_type)) >= len(
            atoms_to_update
        )

        end_atoms_to_update = [atoms_lined_up[0], atoms_lined_up[-1]]
        for atom in end_atoms_to_update:
            assert atom.get_type() != crosslinker_type
            assert atom.get_type() != central_atom_type
            update_atom_type(atom, end_atom_type)
            assert universe.get_atom(atom.get_id()).get_type() == end_atom_type

    # mark one random atom per cluster as "single fixed"
    clusters = universe.get_clusters()
    for cluster in clusters:
        atoms_of_cluster = cluster.get_atoms()
        random.shuffle(atoms_of_cluster)
        for atom in atoms_of_cluster:
            if (
                atom.get_type() == crosslinker_type
                or atom.get_type() == central_atom_type
                or atom.get_type() == end_atom_type
            ):
                continue
            update_atom_type(atom, single_fixed_atom_type)
            break

    masses = {}
    for m in range(max([central_atom_type, end_atom_type, single_fixed_atom_type])):
        masses[m + 1] = 1.0
    universe.set_masses(masses)
    # universe.resample_velocities(0., 1.)
    try:
        print("Temperature: {}".format(universe.compute_temperature()))
    except ValueError:
        pass

    end_atoms = universe.get_atoms_by_type(end_atom_type)
    assert len(end_atoms) >= len(chains)
    central_atoms = universe.get_atoms_by_type(central_atom_type)
    assert len(central_atoms) >= len(chains)
    assert len(central_atoms) == central_atoms_found
    assert n_crosslinks == len(universe.get_atoms_by_type(crosslinker_type))

    file_to_write = input
    # if ("-cem" not in input):
    #     file_to_write = input.replace("structure.out", "structure_cem.out")
    # output new system
    data_writer = DataFileWriter(universe)
    data_writer.config_include_velocities(False)
    data_writer.config_include_angles(True)
    data_writer.config_molecule_idx_for_swap(False)
    data_writer.config_move_into_box(False)
    data_writer.config_attempt_image_reset(True)
    data_writer.write_to_file(file_to_write)
    print('Wrote file "{}"'.format(file_to_write))
