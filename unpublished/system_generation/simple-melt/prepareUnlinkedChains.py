#!/bin/env python

import os
import random
import sys

import numpy as np
import pandas as pd
# https://polymercpp.readthedocs.io/en/latest/
from PolymerCpp.helpers import getCppSAWLC

# input
beadsPerChain = 28
numberOfChains = int(1e4) #int(1e2) # 1e4
beadDistance = 0.96
seed = 8804 + beadsPerChain + numberOfChains
mode_swappable = True # whether to number the molecules per chain or per position in chain

random.seed(seed)
chainsDf = pd.DataFrame(columns=list('xyz'))

originDisplacement = numberOfChains/2
# generate each chain
# NOTE: currently, all chains have the same length
for i in range(numberOfChains):
    # do actual random walking
    chain = getCppSAWLC(beadsPerChain-1, beadDistance, 0.75)
    assert(len(chain) == beadsPerChain)
    chainDf = pd.DataFrame(chain, columns=list('xyz'))
    chainDf[list('xyz')] += (i - originDisplacement) * \
        (1-2*np.array([random.random(), random.random(), random.random()]))
    # append this one chain to all others
    chainsDf = chainsDf.append(chainDf, ignore_index=True)

minima = chainsDf.min(axis=0)  # - 0.2
maxima = chainsDf.max(axis=0)  # + 0.2
# add free space
minima = minima - 3
maxima = maxima + 3  # + 0.5 * lengths
lengths = maxima - minima
# output
print("LAMMPS FENE chain data file: {} chains with {} beads per chain and a momomer distance of {}".format(
    numberOfChains, beadsPerChain, beadDistance))
print("""
          {} atoms
          {} bonds
          {}  angles
          0  dihedrals
          0  impropers

          2  atom types
          1  bond types
          1  angle types
          0  dihedral types
          0  improper types
""".format(len(chainsDf), numberOfChains*(beadsPerChain-1), numberOfChains*(beadsPerChain-2)))
print("""
      {}       {} xlo xhi
      {}       {} ylo yhi
      {}       {} zlo zhi
""".format(minima["x"], maxima["x"], minima["y"], maxima["y"], minima["z"], maxima["z"]))
print("""
Masses

  1  1.0
  2  1.0

Atoms
""")
atomIdx = 1
moleculeIdx = 1
atomType = 1
mIdxCountDirection = "up"
for index, atom in chainsDf.iterrows():
    print("         {}       {}       {}    {}    {}    {}   0   0   0".format(
        atomIdx, moleculeIdx, atomType, atom[0], atom[1], atom[2]))
    atomIdx += 1
    if (mode_swappable):
        if (moleculeIdx == beadsPerChain/2 and mIdxCountDirection == "up"):
            mIdxCountDirection = "none"
        if (moleculeIdx > beadsPerChain/2 and mIdxCountDirection == "up"):
            mIdxCountDirection = "down"
        if (mIdxCountDirection == "up"):
            moleculeIdx += 1
        elif (mIdxCountDirection == "none"):
            mIdxCountDirection = "down"
        else:
            moleculeIdx -= 1
    if ((atomIdx-1) % beadsPerChain == 0):
        if (mode_swappable):
            moleculeIdx = 1
            mIdxCountDirection = "up"
        else:
            moleculeIdx += 1

print("""
Bonds
""")

bondIdx = 1
for i in range(1, beadsPerChain*numberOfChains, 1):
    # skip certain
    if (bondIdx % beadsPerChain != 0):
        print("        {}  1        {}        {}".format(
            bondIdx, bondIdx, bondIdx+1))
    bondIdx += 1

print("""
Angles
""")

angleIdx = 1
for i in range(1, beadsPerChain*numberOfChains, 1):
    # no angle exists between chains
    if ((i % beadsPerChain) != 0 and ((i+1) % beadsPerChain != 0)):
        print("        {}  1        {}        {}        {}".format(
            angleIdx, i, i+1, i+2))
        angleIdx += 1

assert(angleIdx == numberOfChains*(beadsPerChain-2)+1)
