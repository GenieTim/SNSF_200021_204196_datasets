import os
import random
import warnings
from math import sqrt

import matplotlib.pyplot as plt
import numpy as np

# This script is more complicated than necessary,
# as it is originating from the hexagonal one
#
# we first generate the positions of the crosslinkers
# (the ones that actually form the hexagonal lattice)

chainLength = 5**2  # 5 # 5**2
nrOfCellsInDirX = 15  # int(1.5*36) # 4
nrOfCellsInDirY = 15  # 36 # 2
nrOfCellsInDirZ = 15
plotBoundary = True  # False
plotUnitCell = False
includeGhostBonds = False
mc = True
random.seed(7)

# a: chain length, b: box side length, c: half side length
a = chainLength

boxLen = (nrOfCellsInDirX)*a
boxHeight = (nrOfCellsInDirY)*a
boxDepth = (nrOfCellsInDirZ)*a

# setup the crosslinker coordinates
crosslinkerCoordinates = []
bonds = []
currentY = 0
for ycell in range(int(nrOfCellsInDirY)):
    currentX = 0
    for xcell in range(int(nrOfCellsInDirX)):
        currentZ = 0
        for zcell in range(int(nrOfCellsInDirZ)):
            # front face
            crosslinkerCoordinates.extend([
                (currentX, currentY, currentZ)
            ])
            # next it
            currentZ += a
        currentX += a
    currentY += a

# def removeDuplicates(lst):
#     return [t for t in (set(tuple(i) for i in lst))]

# crosslinkerCoordinates = removeDuplicates(crosslinkerCoordinates)


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
    deltaz = abs(coords2[2] - coords1[2])
    if (deltaz > 0.5*boxDepth):
        deltaz -= boxDepth
        isPBC = True
    return np.sqrt(deltax**2 + deltay**2 + deltaz**2), isPBC


# setup the bonds
# they could easily be setup in an algorithm of O(N),
# yet the following is O(N^2). As it only has to run once,
# and is fast enough for our purposes, this comment
# shall serve as a warning if you want to use the algorithm
# for other purposes
availableSites = np.zeros(len(crosslinkerCoordinates)) + 6
maxDist = chainLength if not mc else chainLength*sqrt(chainLength)
for i in range(len(crosslinkerCoordinates)):
    if (availableSites[i] <= 0):
        continue
    potentialBonds = []
    indicesToCompare = list(range(0, len(crosslinkerCoordinates)))
    random.shuffle(indicesToCompare)
    for j in indicesToCompare:
        if (availableSites[j] <= 0):
            continue
        if (availableSites[i] <= 0):
            break
        if (i == j and availableSites[i] < 2):
            continue
        if (not mc and i == j):
            continue
        delta, isPeriodic = computeDistanceBetweenCoords(
            crosslinkerCoordinates[i], crosslinkerCoordinates[j])
        if (delta < 1.005*maxDist):
            sortedRes = sorted([i, j])
            availableSites[i] -= 1
            availableSites[j] -= 1
            sortedRes.append(isPeriodic)
            bonds.append(tuple(sortedRes))

if (np.sum(availableSites) != 0.0):
    warnings.warn("Remaining available sites: {}".format(
        np.sum(availableSites)))
# remove duplicate bonds (every one should be duplicate currently)
bonds = list(set(bonds))

# fig, ax = plt.subplots(projection='3d')
fig = plt.figure()
ax = fig.add_subplot(1, 1,1, projection='3d')
# fig.set_size_inches(boxLen/chainLength*0.875*1.5, boxHeight /
#                     chainLength*0.875*1.5, forward=True)
fig.set_size_inches(7.0, 6.062177826491079, forward=True)

crosslinkerCoordinates = np.array(crosslinkerCoordinates)

# the offset of the box we draw; bonds cannot be split in same lengths
boxOffsetY = -(0.0)*a
boxOffsetX = -(0.25)*a
boxCoordsY = [boxOffsetY, boxHeight+boxOffsetY]
boxCoordsX = [boxOffsetX, boxLen+boxOffsetX]


def computePeriodicAtomNeighbour(atomFrom, atomTo):
    newAtom = atomTo.copy()
    deltax = abs(atomFrom[0] - atomTo[0])
    if (deltax > 0.5*boxLen):
        if (atomTo[0] > 0.5*boxLen):
            newAtom[0] -= boxLen
        else:
            newAtom[0] += boxLen
    deltay = abs(atomFrom[1] - atomTo[1])
    if (deltay > 0.5*boxHeight):
        if (atomTo[1] > 0.5*boxHeight):
            newAtom[1] -= boxHeight
        else:
            newAtom[1] += boxHeight
    deltaz = abs(atomFrom[2] - atomTo[2])
    if (deltaz > 0.5*boxDepth):
        if (atomTo[2] > 0.5*boxDepth):
            newAtom[2] -= boxDepth
        else:
            newAtom[2] += boxDepth
    return newAtom


# plot cross-linker atoms
ax.scatter3D(crosslinkerCoordinates[:, 0],
             crosslinkerCoordinates[:, 1],
             crosslinkerCoordinates[:, 2], s=89, zorder=2)

# plot box
lineThickness = 3
# ax.plot([boxCoordsX[0], boxCoordsX[0]], [boxCoordsY[0], boxCoordsY[1]],  '-',
#         color='tab:orange', linewidth=lineThickness, zorder=5)
# ax.plot([boxCoordsX[0], boxCoordsX[1]], [boxCoordsY[1], boxCoordsY[1]], '-',
#         color='tab:orange', linewidth=lineThickness, zorder=5)
# ax.plot([boxCoordsX[1], boxCoordsX[1]], [boxCoordsY[0], boxCoordsY[1]], '-',
#         color='tab:orange', linewidth=lineThickness, zorder=5)
# ax.plot([boxCoordsX[0], boxCoordsX[1]], [boxCoordsY[0], boxCoordsY[0]], '-',
#         color='tab:orange', linewidth=lineThickness, zorder=5)

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
        ax.plot3D([atomFrom[0], atomTo[0]], [atomFrom[1], atomTo[1]], [atomFrom[2], atomTo[2]], '-',
                  color="black", zorder=1)
    else:
        pass
        continue
        # plot the "cut-off", periodic bonds
        def plotLineHalfDashed(afrom, ato, color="black"):
            x0 = afrom[0]
            x1 = ato[0]
            if (x1 == x0):
                x1 = x0 + 0.5*a
            y0 = afrom[1]
            y1 = ato[1]
            if (y1 == y0):
                y1 = x0 + 0.5*a
            slope = (y0-y1)/(x0-x1)
            intercept = y1 - slope*x1

            def y(x):
                return slope*x + intercept

            def x(y):
                if (slope == 0):
                    return x1
                return (y-intercept)/slope

            def findNearestPoint(origin, target, pointsToCheck):
                nearestPoint = None
                nearestPointD = None
                for point in pointsToCheck:
                    # TODO: limit: only within box spawned by x0, y0, x1, y1
                    if (point[0] > origin[0] and point[0] > target[0]):
                        continue
                    if (point[0] < origin[0] and point[0] < target[0]):
                        continue
                    if (point[1] > origin[1] and point[1] > target[1]):
                        continue
                    if (point[1] < origin[1] and point[1] < target[1]):
                        continue

                    currentPointD = (point[0]-origin[0]
                                      )**2 + (point[1]-origin[1])**2
                    if (nearestPointD is None or currentPointD < nearestPointD):
                        nearestPointD = currentPointD
                        nearestPoint = point

                if (nearestPoint is None):
                    nearestPoint = target

                return nearestPoint

            # find the actually relevant point
            pointOfIntercept = findNearestPoint((x0, y0), (x1, y1), [
                (boxCoordsX[0], y(boxCoordsX[0])),
                (boxCoordsX[1], y(boxCoordsX[1])),
                (x(boxCoordsY[0]), boxCoordsY[0]),
                (x(boxCoordsY[1]), boxCoordsY[1])])

            # solid part
            ax.plot([x0, pointOfIntercept[0]], [
                y0, pointOfIntercept[1]], '-', color=color, zorder=1)
            # dashed part
            ax.plot([pointOfIntercept[0], x1], [
                pointOfIntercept[1], y1], ':', color=color, zorder=1)

        # need to plot twice the same bond, half
        atomNewTo = computePeriodicAtomNeighbour(atomFrom, atomTo)
        atomNewFrom = computePeriodicAtomNeighbour(atomTo, atomFrom)

        plotLineHalfDashed(atomFrom, atomNewTo)
        plotLineHalfDashed(atomTo, atomNewFrom)

for i in range(len(crosslinkerCoordinates)):
    ax.text(crosslinkerCoordinates[i][0], crosslinkerCoordinates[i][1], crosslinkerCoordinates[i][2], str(i), 'x')

if (not plotBoundary):
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    ax.set_frame_on(False)

# plot unit cell
if (plotUnitCell):
    unitCellLineThickness = 2  # lineThickness/2
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
            "figure{}-3d-primitive-cubic-lattice-{}{}-{}x{}x{}.png".format("-mc" if mc else "", "no-boundary-" if not plotBoundary else "", chainLength, nrOfCellsInDirX, nrOfCellsInDirY, nrOfCellsInDirZ)), bbox_inches="tight", dpi=300)

# verify the functionality of each atom
for i in range(len(crosslinkerCoordinates)):
    if (crosslinkerFunctionality[i] != 6):
        warnings.warn("Found crosslinker {} with functionality {} != 6".format(
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
    deltaZ = atomTo[2] - atomFrom[2]
    currentX = atomFrom[0]
    currentY = atomFrom[1]
    currentZ = atomFrom[2]
    for newAtomI in range(chainLength-1):
        currentX += deltaX/chainLength
        currentY += deltaY/chainLength
        currentZ += deltaZ/chainLength
        otherAtomCoordinates.append([
            currentX, currentY, currentZ
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
print("LAMMPS pc lattice data file: nx = {}, ny = {}, chain len = {}".format(
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
""".format(0, boxLen, 0, boxHeight, 0, boxDepth))
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
        atomIdx, moleculeIdx, 2, crosslinkerCoordinates[i][0], crosslinkerCoordinates[i][1], crosslinkerCoordinates[i][2]))
    atomIdx += 1

for i in range(len(otherAtomCoordinates)):
    print("         {}       {}       {}    {}    {}    {}   0   0   0".format(
        atomIdx, moleculeIdx, 1, otherAtomCoordinates[i][0], otherAtomCoordinates[i][1], otherAtomCoordinates[i][2]))
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


if (includeGhostBonds):
    print("""
    Bond Coeffs

    1 harmonic 4.1419333334 0
    """)
    print("2 zero")
