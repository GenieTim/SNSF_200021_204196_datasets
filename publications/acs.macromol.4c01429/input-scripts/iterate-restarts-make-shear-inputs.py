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
    # ["shear-everystep_from_restart", "*.restart-*.out"],
    # ["shear-everystep", "t-324100000-t-36000000-crosslinked_M10000_N39_p_0.98.out-equilibration_everystep_cont.structure.out"]
]

output_folder = "/cluster/scratch/betim/doctorate-hobbies-3/md-experiment-simulations/output"

for inputStructure in inputStructures:
    baseDir = os.path.dirname(__file__)
    fileToCopy = os.path.join(baseDir, "input/{}.in".format(inputStructure[0]))

    restartFiles = glob.glob(os.path.join(
        baseDir, "output/{}".format(inputStructure[1])))
    actualRestartFiles = {}
    for file in restartFiles:
        structureBasename = os.path.basename(file).replace(
            ".restart-1", "").replace(".restart-2", "")
        if (structureBasename not in actualRestartFiles):
            actualRestartFiles[structureBasename] = []
        actualRestartFiles[structureBasename].append(file)

    restartFiles = []
    for key, files in actualRestartFiles.items():
        if (len(files) == 1):
            restartFiles.append(files[0])
        else:
            sortedFiles = sorted(
                files, key=lambda x: os.path.getmtime(x), reverse=True)
            restartFiles.append(files[0])  # the newest one

    if (len(restartFiles) == 0):
        warnings.warn(
            "No restart files found for {}".format(inputStructure[1]))
        continue

    with open(fileToCopy, 'r') as f:
        content = f.read()
        for restartFile in restartFiles:
            for dir in ["xy", "xz", "yz"]:
                for gamma in [-0.1, 0.1]:
                    structureBasename = os.path.basename(restartFile).replace(
                        ".restart-1", "").replace(
                            ".restart-2", "").replace(
                                ".out-equilibration_everystep.out", ".out").replace(
                                    ".out-equilibration_everystep_cont.out", ".out")
                    N = extractNFromFilename(structureBasename)
                    inputBasename = os.path.basename(fileToCopy)
                    outputBasename = "{}-{}-{}-{}".format(dir, gamma,
                                                          structureBasename, inputBasename)
                    to_replace = {
                        "&structure-basename": structureBasename,
                        "&output-folder": output_folder,
                        "&structure-file": "./structure/" + structureBasename,
                        "&output-basename": outputBasename.removesuffix(".in"),
                        "&deform-dir": dir,
                        "&deform-rate": str(np.sign(gamma)*0.001),
                        "&deform-timestep": "0.001",
                        "&deform-steps": str(int(abs(gamma)/(0.001*0.001))),
                        "&restart-file": os.path.join(output_folder, os.path.basename(restartFile))
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
                    print("Wrote {}".format(target_file))
