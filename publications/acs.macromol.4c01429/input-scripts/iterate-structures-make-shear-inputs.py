#!/usr/bin/env python
import glob
import os
import warnings

import numpy as np
from pylimer_doctorate_utils.fileConfigResolver import extractNFromFilename
from pylimer_doctorate_utils.myre import preg_match_group

inputStructures = [
    # ["shear", "*.out"],
    # ["shear-short", "*N19_p_0.98*.out"],
    ["shear-everystep", "*t-237100000-t-123600000-t-57700000-crosslinked_M10000_N79_p_0.95*"],
    # ["shear-everystep", "t-324100000-t-36000000-crosslinked_M10000_N39_p_0.98.out-equilibration_everystep_cont.structure.out"]
]

for inputStructure in inputStructures:
    baseDir = os.path.dirname(__file__)
    fileToCopy = os.path.join(baseDir, "input/{}.in".format(inputStructure[0]))

    structure_files = glob.glob(os.path.join(
        baseDir, "structure/{}".format(inputStructure[1])))
    structure_files = [
        file for file in structure_files if not "_type_per_molecule" in file]

    if (len(structure_files) == 0):
        warnings.warn(
            "No structure files found for {}".format(inputStructure[1]))
        continue

    with open(fileToCopy, 'r') as f:
        content = f.read()
        for structure in structure_files:
            for dir in ["xy", "xz", "yz"]:
                for gamma in [-0.1, 0.1]:
                    structureBasename = os.path.basename(structure)
                    N = extractNFromFilename(structureBasename)
                    inputBasename = os.path.basename(fileToCopy)
                    outputBasename = "{}-{}-{}-{}".format(dir, gamma,
                                                          structureBasename, inputBasename)
                    to_replace = {
                        "&structure-basename": structureBasename,
                        "&output-folder": "/cluster/scratch/betim/doctorate-hobbies-3/md-experiment-simulations/output",
                        "&structure-file": "./structure/" + structureBasename,
                        "&output-basename": outputBasename.removesuffix(".in"),
                        "&deform-dir": dir,
                        "&deform-rate": str(np.sign(gamma)*0.001),
                        "&deform-timestep": "0.001",
                        "&deform-steps": str(int(abs(gamma)/(0.001*0.001)))
                    }
                    replacedContent = content
                    for key, val in to_replace.items():
                        replacedContent = replacedContent.replace(key, val)
                    if ("&" in replacedContent):
                        lines = replacedContent.split("\n")
                        for line in lines:
                            if ("&" in line):
                                print("Not-replaced in line: '{}'".format(line))
                    # make sure we replaced everything we need
                    assert("&" not in replacedContent)
                    target_file = os.path.join(
                        baseDir, "input-generated", "{}".format(outputBasename))
                    with open(target_file, 'w') as f:
                        f.write(replacedContent)
