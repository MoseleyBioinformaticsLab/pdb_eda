"""
Utilities (pdb_eda.utils)
-------------------------

Contains low-level functions used in pdb_eda.ccp4 and pdb_eda.densityAnalysis.
Show not be used, but is a fallback if cutils.pyx cannot be cythonized.
The cythonized version provide a 3- to 4-fold improvement in execution performance.
"""
def testOverlap(selfBlob, otherBlob):
    """Check if two blobs overlaps or right next to each other.

    :param selfBlob: A :class:`pdb_eda.ccp4.DensityBlob` object.
    :param otherBlob: A :class:`pdb_eda.ccp4.DensityBlob` object.
    :return: :py:obj:`True` or :py:obj:`False`.
    """
    return True if any(-1 <= x[0] - y[0] <= 1 and -1 <= x[1] - y[1] <= 1 and -1 <= x[2] - y[2] <= 1 for x in selfBlob.crsList for y in otherBlob.crsList) else False


def sumOfAbs(array, cutoff):
    """Return sum of absolute values above a cutoff.

    :param iterable array:
    :param float cutoff:
    :return: value
    :rtype: float
    """
    return sum(abs(value) for value in array if abs(value) > cutoff)

import numpy as np
import scipy.spatial
dcutoff = np.sqrt(3)  ## the points are considered to be adjacent if one is in the one layer outer box with the other one in the center
def createCrsLists(crsList):
    """Calculates a list of crsLists from a given crsList.
    This is a preparation step for creating blobs.

    :param crsList: a crs list.
    :return: crsLists is a list of crsList.
    :rtype: :py:obj:`list`
    """
    crsArray = np.matrix(crsList)
    distances = scipy.spatial.distance.cdist(crsArray, crsArray)

    crsLists = []
    usedIdx = set()
    for startingIndex in range(len(crsList)):
        if not startingIndex in usedIdx:
            currCluster = set([startingIndex])
            newCluster = {index for index, distance in enumerate(distances[startingIndex]) if index not in currCluster and distance <= dcutoff}
            currCluster.update(newCluster)
            while len(newCluster):
                newCluster = {index for oldIndex in newCluster for index, distance in enumerate(distances[oldIndex]) if index not in currCluster and distance <= dcutoff}
                currCluster.update(newCluster)

            usedIdx.update(currCluster)
            crsLists.append([crsList[index] for index in currCluster])
    return crsLists

import itertools
def createSymmetryAtoms(atomList, rotationMats, orthoMat, xs, ys, zs):
    """Creates and returns a list of all symmetry atoms.

    :param :py:obj:`list` atomList:
    :param list rotationMats:
    :param list orthoMat:
    :param list xs:
    :param list ys:
    :param list zs:
    :return: allAtoms
    :rtype: :py:obj:`list`
    """
    allAtoms = []
    for symmetry in itertools.product([-1, 0, 1],[-1, 0, 1],[-1, 0, 1],range(len(rotationMats))):
        if symmetry == (0,0,0,0):
            allAtoms.extend([SymAtom(atom, atom.coord, symmetry) for atom in atomList])
        else:
            rMat = rotationMats[symmetry[3]]
            otMat = np.dot(orthoMat, symmetry[0:3])
            coordList = [np.dot(rMat[:, 0:3], atom.coord) + rMat[:, 3] + otMat for atom in atomList]
            allAtoms.extend([SymAtom(atom, coord, symmetry) for atom,coord in zip(atomList,coordList)
                if xs[0] - 5 <= coord[0] <= xs[-1] + 5 and ys[0] - 5 <= coord[1] <= ys[-1] + 5 and zs[0] - 5 <= coord[2] <= zs[-1] + 5])

    return allAtoms

class SymAtom:
    """A wrapper class to the `BioPDB.atom` class, delegating all BioPDB atom class methods and data members except having its own symmetry and coordination."""

    def __init__(self, atom, coord, symmetry):
        """
        `pdb_eda.densityAnalysis.symAtom` initializer.

        :param `BioPDB.atom` atom: atom object.
        :param :py:obj:`list` coord: x,y,z coordinates.
        :param :py:obj:`list` symmetry: i,j,k,r symmetry
        """
        self.atom = atom
        self.coord = coord
        self.symmetry = symmetry

    def __getattr__(self, attr):
        return getattr(self.atom, attr)

def getPointDensityFromCrs(densityMatrix, crsCoord):
    """Get the density of a point.

    :param crsCoord: crs coordinates.
    :type crsCoord: A :py:obj:`list` of :py:obj:`int`
    :return: density
    :rtype: float
    """
    crsCoord = list(crsCoord)
    header = densityMatrix.header
    for ind in range(3):
        if crsCoord[ind] < 0 or crsCoord[ind] >= header.ncrs[ind]:
            crsCoord[ind] -= int(np.floor(crsCoord[ind] / header.crsInterval[ind]) * header.crsInterval[ind])

        if header.ncrs[ind] <= crsCoord[ind] < header.crsInterval[ind]: # think this should include "or crsCoord[ind] < 0"
            return 0

    return densityMatrix.density[crsCoord[2], crsCoord[1], crsCoord[0]]

def createFullCrsList(densityMatrix, cutoff):
    """Returns full crs list for the density matrix.

    :param densityMatrix:
    :param float cutoff:
    :return: crsList
    :rtype: :py:obj:`list`
    """
    ## only explore the non-repeating part (<= # xyz intervals) of the density map for blobs
    ncrs = densityMatrix.header.uniqueNcrs
    if cutoff > 0:
        return [ crs for crs in itertools.product(range(ncrs[0]),range(ncrs[1]),range(ncrs[2])) if getPointDensityFromCrs(densityMatrix, crs) >= cutoff ]
    elif cutoff < 0:
        return [ crs for crs in itertools.product(range(ncrs[0]),range(ncrs[1]),range(ncrs[2])) if getPointDensityFromCrs(densityMatrix, crs) <= cutoff ]
    else:
        return None
