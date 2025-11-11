#!/usr/bin/env python
import glob
import os
import random
import warnings

from pylimer_doctorate_utils.fileConfigResolver import extractNFromFilename
from pylimer_doctorate_utils.myre import preg_match_group
from pylimer_tools.io.read_lammps_output_file import read_data_file

input_structures = [
    # ["equilibration_and_more/equilibration_autocorrelation", "*.out"],
    # ["equilibration_and_more/equilibration", "crosslinked/crosslinked_p_1.0_melt_100_a_125_50_xlinks_v_*.out"],
    # ["equilibration_and_more/equilibration_do_crosslink", "melt_100_a_*_100_xlinks_v_*.V-fixed.structure.out"],
    # ["equilibration_and_more/equilibration_long_do_crosslink", "*xlinks_mono*.structure.out"],
    # ["equilibration_and_more/equilibration_long_do_crosslink_w_measurement", "pre-crosslinking/*xlinks_mono*"],
    # ["equilibration_and_more/equilibration_very_slow", "*diamond-lattice*"]
    # ["equilibration", "polydisperse/*p_choice2_melt_100_a_125_*.V-fixed*"],
    # ["equilibration", "melt_100_a_640_v*"],
    # ["equilibration_long_do_crosslink", "melt_100_a_638_50_xlinks_v*"]
    # ["sampling/sampling_gmid_ree_incl_long_no_angle", "*-cem.structure.out"]
    # ["equilibration_and_more/equilibration_sampling_gmid_ree_incl_long_no_angle",
    #     "melt_100_a_3_v*-cem.structure.out"]
    # ["equilibration_and_more/equilibration_sampling_incl_long", "crosslinked/crosslinked_p_1_1_melt_100_a_*_50*"]
    # ["equilibration_and_more/equilibration_sampling_incl_long_xlong_no_angle",
    #     "crosslinked/crosslinked_p_1_0.5_melt_100_a_18_100_xlinks_v_*.structure.out"]
    ["equilibration_and_more/equilibration_do_crosslink", "melt_10000_a_125_5000_xlinks_v_*.structure.out"]
    # [
    #     "equilibration_and_more/equilibration_sampling_incl_long_xlong_no_angle",
    #     "crosslinked/crosslinked_p_1*melt_10000_a_125_8500_xlinks_v_*.structure.out",
    # ]
    # ["equilibration_do_crosslink_polydisperse", "*100_a_125_85_*"]
    # ["equilibration_sampling_gmid_ree_incl_long_no_angle",
    #     "*-cem.structure.out"]
]

output_folder = (
    "/cluster/scratch/betim/doctorate-hobbies-4/md-experiment-simulations/output"
)

for input_structure in input_structures:
    base_dir = os.path.dirname(__file__)
    file_to_copy = os.path.join(base_dir, "input/{}.in".format(input_structure[0]))

    structure_files = glob.glob(
        os.path.join(base_dir, "structure/{}".format(input_structure[1]))
    )
    structure_files = [
        file for file in structure_files if "_type_per_molecule" not in file
    ]

    if len(structure_files) == 0:
        warnings.warn("No structure files found for {}".format(input_structure[1]))
        continue

    with open(file_to_copy, "r") as f:
        content = f.read()
        for structure_file in structure_files:
            v = int(preg_match_group(structure_file, r"v_([0-9]*)[\.\-]"))  # V-fixed
            # if v not in [42, 35, 28, 16, 31, 34, 4, 49, 6, 8, 25, 26, 15, 23, 44, 11, 12, 50, 39, 5, 9, 10, 17, 33, 38]:
            #     continue
            # if (v < 26):
            #     continue
            structure_basename = os.path.basename(structure_file)
            structure_dirname = os.path.dirname(structure_file)
            structure_subdir = structure_dirname.split("/structure", 1)[1]
            N = extractNFromFilename(structure_basename)
            input_basename = os.path.basename(file_to_copy)
            output_basename = "{}-{}".format(structure_basename, input_basename)
            replaced_content = content

            to_replace = {
                "&structure-basename": structure_basename,
                "&structure-file": "./structure"
                + (structure_subdir + "/" if structure_subdir != "" else "/")
                + structure_basename,
                "&output-basename": output_basename.removesuffix(".in"),
                "&output-folder": output_folder,
                "&rand-seed": str(random.randint(0, 100000)),
            }

            for key, value in to_replace.items():
                replaced_content = replaced_content.replace(key, value)

            if "&nchains" in replaced_content:
                Nchains = preg_match_group(structure_basename, r"melt_([0-9]*)_a_")
                replaced_content = replaced_content.replace("&nchains", Nchains)
            if "&nmonochains" in replaced_content:
                Nchains = preg_match_group(structure_basename, r"mono_([0-9]*)_a_")
                if Nchains is None:
                    Nchains = "0"
                replaced_content = replaced_content.replace("&nmonochains", Nchains)

            if "&numBeadsType" in replaced_content:
                universe = read_data_file(structure_file)
                atom_type_count = universe.count_atom_types()
                for type, count in atom_type_count.items():
                    replaced_content = replaced_content.replace(
                        "&numBeadsType" + str(type), str(count)
                    )

            if "&" in replaced_content:
                lines = replaced_content.split("\n")
                for line in lines:
                    if "&" in line:
                        print("Not-replaced in line: '{}'".format(line))
            # make sure we replaced everything we need
            assert "&" not in replaced_content
            target_file = os.path.join(
                base_dir, "input-generated", "{}".format(output_basename)
            )
            with open(target_file, "w") as f:
                f.write(replaced_content)
            print("Wrote {}".format(target_file))
