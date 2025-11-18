from math import sqrt
import os
import warnings
import random

import matplotlib.pyplot as plt
import numpy as np

# we first generate the positions of the crosslinkers
# (the ones that actually form the hexagonal lattice)

chainLength = 5**2 # 5 # 5**2
nrOfCellsInDirX = 4
nrOfCellsInDirY = 4
plotBoundary = False
includeGhostBonds = False
random.seed(5) # 3

# a: side length, b: longer tilted side length, c: shorter tilted side length
a = chainLength
b = np.sqrt(3)/2 * a
c = 0.5*a

# make sure the periodicity can be achieved
if (nrOfCellsInDirX % 2 != 0):
    nrOfCellsInDirX += 1
if (nrOfCellsInDirY % 2 != 0):
    nrOfCellsInDirY += 1

boxLen = (nrOfCellsInDirX*2)*b
boxHeight = (nrOfCellsInDirY)*(a+c)

# setup the crosslinker coordinates
crosslinkerCoordinates = []
bonds = []
currentY = a
for ycell in range(int(nrOfCellsInDirY/2)):
    currentX = 0.5*b
    for xcell in range(int(nrOfCellsInDirX)):
        crosslinkerCoordinates.append([
            currentX, currentY
        ])
        crosslinkerCoordinates.append([
            currentX+b, currentY-c
        ])
        crosslinkerCoordinates.append([
            currentX, currentY+a
        ])
        crosslinkerCoordinates.append([
            currentX+b, currentY+a+c
        ])
        currentX += 2*b
    currentY += 2*(c + a)


def computeDistanceBetweenCoords(coords1, coords2):
    isPBC = False
    deltax = abs(coords2[0] - coords1[0])
    if (deltax > 0.5*boxLen):
        deltax -= boxLen
        isPBC = True
    deltay = abs(coords2[1] - coords1[1])
    if (deltay > 0.5*boxHeight):
        deltay -= boxHeight
        isPBC = True
    return np.sqrt(deltax**2 + deltay**2), isPBC


def computePeriodicAtomNeighbour(atomFrom, atomTo):
    newAtom = atomTo.copy()
    deltax = abs(atomFrom[0] - atomTo[0])
    if (deltax > 0.5*boxLen):
        if (atomTo[0] > 0.5*boxLen+0.5*b):
            newAtom[0] -= boxLen
        else:
            newAtom[0] += boxLen
    deltay = abs(atomFrom[1] - atomTo[1])
    if (deltay > 0.5*boxHeight):
        if (atomTo[1] > 0.5*boxHeight+a):
            newAtom[1] -= boxHeight
        else:
            newAtom[1] += boxHeight
    return newAtom


# setup the bonds
# they could easily be setup in an algorithm of O(N),
# yet the following is O(N^2). As it only has to run once,
# and is fast enough for our purposes, this comment
# shall serve as a warning if you want to use the algorithm
# for other purposes
availableSites = np.zeros(len(crosslinkerCoordinates)) + 3
for i in range(len(crosslinkerCoordinates)):
    if (availableSites[i] <= 0):
        continue
    potentialBonds = []
    indicesToCompare = list(range(0, len(crosslinkerCoordinates)))
    random.shuffle(indicesToCompare)
    for j in indicesToCompare:
        if (availableSites[i] <= 0):
            break
        if (availableSites[j] <= 0):
            continue
        if (i == j and availableSites[i] < 2):
            continue
        delta, isPeriodic = computeDistanceBetweenCoords(
            crosslinkerCoordinates[i], crosslinkerCoordinates[j])
        if (delta < 1.005*chainLength*sqrt(chainLength)):
            sortedRes = sorted([i, j])
            availableSites[i] -= 1
            availableSites[j] -= 1
            sortedRes.append(isPeriodic)
            bonds.append(tuple(sortedRes))

# remove duplicate bonds (every one should be duplicate currently)
bonds = list(set(bonds))

fig, ax = plt.subplots()
fig.set_size_inches(boxLen/chainLength*0.875, boxHeight /
                    chainLength*0.875, forward=True)

crosslinkerCoordinates = np.array(crosslinkerCoordinates)

# count functionality for verification purposes
crosslinkerFunctionality = {}
# plot bonds (and count functionality)
for i in range(len(bonds)):
    bond = bonds[i]
    if (bond[0] not in crosslinkerFunctionality):
        crosslinkerFunctionality[bond[0]] = 0
    crosslinkerFunctionality[bond[0]] += 1

    if (bond[1] not in crosslinkerFunctionality):
        crosslinkerFunctionality[bond[1]] = 0
    crosslinkerFunctionality[bond[1]] += 1

    atomFrom = crosslinkerCoordinates[bond[0]]
    atomTo = crosslinkerCoordinates[bond[1]]
    if (not bond[2]):
        ax.plot([atomFrom[0], atomTo[0]], [atomFrom[1], atomTo[1]], '-',
                color="black", zorder=1)
    else:
        # plot the "cut-off", periodic bonds
        def plotLineHalfDashed(afrom, ato, color="black"):
            x0 = afrom[0]
            x1 = ato[0]
            y0 = afrom[1]
            y1 = ato[1]
            ax.plot([x0, x0+(x1-x0)/2], [
                y0, y0+(y1-y0)/2], '-', color=color, zorder=1)
            # dashed part
            ax.plot([x0+(x1-x0)/2, x1], [
                y0+(y1-y0)/2, y1], ':', color=color, zorder=1)

        # need to plot twice the same bond, half
        atomNewTo = computePeriodicAtomNeighbour(atomFrom, atomTo)
        atomNewFrom = computePeriodicAtomNeighbour(atomTo, atomFrom)

        plotLineHalfDashed(atomFrom, atomNewTo)
        plotLineHalfDashed(atomTo, atomNewFrom)

# plot atoms
ax.scatter(crosslinkerCoordinates[:, 0],
           crosslinkerCoordinates[:, 1], s=89, zorder=2)

# for i in range(len(crosslinkerCoordinates)):
#     ax.annotate(str(i), crosslinkerCoordinates[i])

# plot box
lineThickness = 3
ax.plot([0, 0], [0, boxHeight],  '-', color='tab:orange', linewidth=lineThickness)
ax.plot([0, boxLen], [boxHeight, boxHeight], '-', color='tab:orange', linewidth=lineThickness)
ax.plot([boxLen, boxLen], [0, boxHeight], '-', color='tab:orange', linewidth=lineThickness)
ax.plot([0, boxLen], [0, 0], '-', color='tab:orange', linewidth=lineThickness)
if (not plotBoundary):
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    ax.set_frame_on(False)

# plot unit cell
unitCellLineThickness = 2#lineThickness/2
smallAdjust = 0.0  # 0.25
ax.plot([0+smallAdjust, 0+smallAdjust],
        [0+smallAdjust, 3*a+smallAdjust], '--', color='tab:green', linewidth=unitCellLineThickness)
ax.plot([0+smallAdjust, 2*b+smallAdjust],
        [0+smallAdjust, 0+smallAdjust], '--', color='tab:green', linewidth=unitCellLineThickness)
ax.plot([2*b+smallAdjust, 2*b+smallAdjust],
        [0+smallAdjust, 3*a+smallAdjust], '--', color='tab:green', linewidth=unitCellLineThickness)
ax.plot([0+smallAdjust, 2*b+smallAdjust],
        [3*a+smallAdjust, 3*a+smallAdjust], '--', color='tab:green', linewidth=unitCellLineThickness)

# ax.set(title="len: {}, nX: {}, nY: {}".format(
#     chainLength, nrOfCellsInDirX, nrOfCellsInDirY))
fig.savefig(os.path.join(os.path.dirname(__file__),
            "figure-2d-hexagonal-mc-lattice-{}{}-{}x{}.png".format("no-boundary-" if not plotBoundary else "", chainLength, nrOfCellsInDirX, nrOfCellsInDirY)))

# verify the functionality of each atom
for i in range(len(crosslinkerCoordinates)):
    if (crosslinkerFunctionality[i] != 3):
        warnings.warn("Found crosslinker {} with functionality {} != 3".format(
            i, crosslinkerFunctionality[i]))

# add "in-between" atoms on every bond
actualBonds = []
otherAtomCoordinates = []
basisId = len(crosslinkerCoordinates)
for i in range(len(bonds)):
    bond = bonds[i]
    atomFrom = crosslinkerCoordinates[bond[0]]
    atomTo = crosslinkerCoordinates[bond[1]]
    if (bond[2]):
        # periodic case: we cheat a bit. Let the atoms escape the box, they will be wrapped by LAMMPS
        atomTo = computePeriodicAtomNeighbour(atomFrom, atomTo)
    # non-periodic case: "easy"
    deltaX = atomTo[0] - atomFrom[0]
    deltaY = atomTo[1] - atomFrom[1]
    currentX = atomFrom[0]
    currentY = atomFrom[1]
    for newAtomI in range(chainLength-1):
        currentX += deltaX/chainLength
        currentY += deltaY/chainLength
        otherAtomCoordinates.append([
            currentX, currentY
        ])
        thisAtomId = basisId+len(otherAtomCoordinates)-1
        if (newAtomI == 0):
            # opening one
            actualBonds.append([bond[0], thisAtomId])
        else:
            actualBonds.append([thisAtomId-1, thisAtomId])
    # and the closing one
    actualBonds.append([bond[1], basisId+len(otherAtomCoordinates)-1])

# finally, write the LAMMPS data file
print("LAMMPS hexagonal lattice data file: nx = {}, ny = {}, chain len = {}".format(
    nrOfCellsInDirX, nrOfCellsInDirY, chainLength))
print("""
          {} atoms
          {} bonds
          0  angles
          0  dihedrals
          0  impropers

          2  atom types
          {}  bond types
          0  angle types
          0  dihedral types
          0  improper types
""".format(len(otherAtomCoordinates)+len(crosslinkerCoordinates), len(actualBonds)+(len(bonds) if includeGhostBonds else 0), 2 if includeGhostBonds else 1))
print("""
      {}       {} xlo xhi
      {}       {} ylo yhi
      {}       {} zlo zhi
""".format(0, boxLen, 0, boxHeight, -0.5, 0.5))
print("""
Masses

  1  1.0
  2  1.0

Atoms
""")
moleculeIdx = 1
atomIdx = 1
for i in range(len(crosslinkerCoordinates)):
    print("         {}       {}       {}    {}    {}    {}   0   0   0".format(
        atomIdx, moleculeIdx, 2, crosslinkerCoordinates[i][0], crosslinkerCoordinates[i][1], 0))
    atomIdx += 1

for i in range(len(otherAtomCoordinates)):
    print("         {}       {}       {}    {}    {}    {}   0   0   0".format(
        atomIdx, moleculeIdx, 1, otherAtomCoordinates[i][0], otherAtomCoordinates[i][1], 0))
    atomIdx += 1

print("""
Bonds
""")

bondIdx = 1
for i in range(len(actualBonds)):
    print("        {}  1        {}        {}".format(
        i, actualBonds[i][0]+1, actualBonds[i][1]+1))
    bondIdx += 1

# also print the bonds of the crosslinkers as zero bonds
if (includeGhostBonds):
    for i in range(len(bonds)):
        print("        {}  2        {}        {}".format(
            i, bonds[i][0]+1, bonds[i][1]+1))
        bondIdx += 1


print("""
Bond Coeffs

1 harmonic 4.1419333334 0
""")
if (includeGhostBonds):
    print("2 zero")


