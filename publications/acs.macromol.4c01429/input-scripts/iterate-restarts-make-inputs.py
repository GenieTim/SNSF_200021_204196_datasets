#!/usr/bin/env python
import glob
import os
import random
import warnings

from pylimer_doctorate_utils.fileConfigResolver import (
    extractNFromFilename, findStructureFileForOutput)
from pylimer_doctorate_utils.myre import preg_match_group
from pylimer_tools.io.read_lammps_output_file import read_data_file

input_structures = [
    # ["equilibration/equilibration_autocorrelation", "*.out"],
    # ["equilibration/equilibration_everystep_cont", "*320*.restart-*.out"],
    # ["crosslinking/do_crosslink", "*equilibration_do_crosslink*.restart-*.out"],
    # ["equilibration/equilibration", "*.out"],
    # ["sampling/sampling_certain_xlong", "melt_100_a_320*sampling_xlong.restart-*.out"],
    # ["sampling/sampling_mlong_T_1.25", "melt_100_a_160_v_*.structure.out-equilibration.out-sampling.restart-*.out"]
    # ["sampling/sampling_xxlong", "crosslinked_p_0.98_melt_100_a_125_50_xlinks_*sampling_xlong.restart*"]
    # ["crosslinking/continue_do_crosslink_polydisperse",
    #     "melt_100_a_125_85_xlinks_v_*.V-fixed.structure.out-equilibration_do_crosslink_polydisperse.restart*"]
    # ["sampling/sampling_long", "melt_100_a_80_v_*sampling.restart*"]
    # ["sampling/sampling_long_split", "crosslinked_p_0.98_melt_100_a_18_50_xlinks_v_*sampling.restart*"]
    # ["sampling/sampling_xlong", "melt_100_a_125_v*-sampling_long.restart-*.out"]
    ["crosslinking/finish_crosslinking", "*_10000_a_125_8500_xlinks_*.restart*"]
    # ["sampling/sampling_xlong_no_angle", "crosslinked_p_1_1_*.restart*"]
    # ["sampling/sampling_xxlong_no_angle", "crosslinked_p_1*_melt_100_a_18_*xlong*no_angle*.restart*"]
    # ["sampling/sampling_gmid_ree_incl_long_no_angle",
    # "melt_100_a_160_v_*.V-fixed-cem.structure.out-equilibration_sampling_gmid_ree_incl_long_no_angle.restart*.out"]
    # ["crosslinking/do_crosslink_polydisperse",
    #     "melt_100_a_638_50_xlinks_v_*.V-fixed.structure.out-equilibration_long_do_crosslink.restart-*.out"]
    # ["sampling/sampling_gmid_ree_incl_long_no_angle",
    #     "melt_100_a_640_v_*.V-fixed.structure.out-equilibration.restart-*.out"]
    # ["equilibration_and_more/restart_to_data",
    #     "melt_100_a_640_v_*.V-fixed.structure.out-equilibration.restart-*.out"]
    # ["sampling/sampling_gmid_incl_long_xlong_no_angle", "*.out-equilibration_sampling_gmid_incl_long_xlong_no_angle.restart*"]
]

output_folder = "/cluster/scratch/betim/doctorate-hobbies-4/md-experiment-simulations/output"

random.seed(8003)


base_dir = os.path.dirname(__file__)
all_structure_files = glob.glob(os.path.join(
    base_dir, "structure/**/*.structure.out"), recursive=True)

for input_structure in input_structures:
    file_to_copy = os.path.join(
        base_dir, "input/{}.in".format(input_structure[0]))

    restart_files = glob.glob(os.path.join(
        base_dir, "output/{}".format(input_structure[1])))
    actual_restart_files = {}
    for file in restart_files:
        file_key = os.path.basename(file).replace(
            ".restart-1", "").replace(".restart-2", "")
        if (file_key not in actual_restart_files):
            actual_restart_files[file_key] = []
        actual_restart_files[file_key].append(file)

    restart_files = []
    for key, files in actual_restart_files.items():
        if (len(files) == 1):
            restart_files.append(files[0])
        else:
            sorted_files = sorted(
                files, key=lambda x: os.path.getmtime(x), reverse=True)
            restart_files.append(files[0])  # the newest one

    if (len(restart_files) == 0):
        warnings.warn(
            "No restart files found for {}".format(input_structure[1]))
        continue

    with open(file_to_copy, 'r') as f:
        content = f.read()
        for restart_file in restart_files:
            restart_basename = os.path.basename(restart_file)
            structure_file = findStructureFileForOutput(
                restart_file, all_structure_files)
            structure_basename = os.path.basename(structure_file)
            # .replace(".out-equilibration_do_crosslink_polydisperse", "")
            N = extractNFromFilename(structure_basename)
            input_basename = os.path.basename(file_to_copy)
            output_basename = "{}-{}".format(structure_basename,
                                             input_basename)
            replaced_content = content.replace(
                "&restart-file", os.path.join(output_folder, os.path.basename(restart_file)))

            to_replace = {
                "&structure-file": "./structure/" + structure_basename,
                "&output-basename": output_basename.removesuffix(".in"),
                "&previous-basename": restart_basename.replace(
                    ".restart-1.out", "").replace(".restart-2.out", ""),
                "&output-folder": output_folder,
                "&nAtomsPerChain": str(N),
                "&rand-seed": str(random.randint(0, 100000))
            }

            for key, value in to_replace.items():
                replaced_content = replaced_content.replace(key, value)

            if ("&nchains" in replaced_content):
                Nchains = preg_match_group(
                    structure_basename, r"melt_([0-9]*)_a_")
                replaced_content = replaced_content.replace(
                    "&nchains", Nchains)
            if ("&nxlinks" in replaced_content):
                Nxlinks = preg_match_group(
                    structure_basename, r"_([0-9]*)_xlinks")
                replaced_content = replaced_content.replace(
                    "&nxlinks", Nxlinks)
            if ("&nmonochains" in replaced_content):
                Nchains = preg_match_group(
                    structure_basename, r"mono_([0-9]*)_a_")
                if (Nchains is None):
                    Nchains = "0"
                replaced_content = replaced_content.replace(
                    "&nmonochains", Nchains)

            if ("&numBeadsType" in replaced_content):
                universe = read_data_file(structure_file)
                atom_type_count = universe.countAtomTypes()
                for type, count in atom_type_count.items():
                    replaced_content = replaced_content.replace(
                        "&numBeadsType" + str(type), str(count))

            if ("&" in replaced_content):
                lines = replaced_content.split("\n")
                for line in lines:
                    if ("&" in line):
                        print("Not-replaced in line: '{}'".format(line))
            # make sure we replaced everything we need
            assert ("&" not in replaced_content)
            target_file = os.path.join(
                base_dir, "input-generated", "{}".format(output_basename))
            with open(target_file, 'w') as f:
                f.write(replaced_content)
            print("Wrote {}".format(target_file))
