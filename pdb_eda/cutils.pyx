# Compile using:
# $ python3 cutils_setup.py build_ext --inplace
# $ cp build/lib.linux-x86_64-3.6/pdb_eda/cutils.cpython-36m-x86_64-linux-gnu.so cutils.so

def testOverlap(selfBlob, otherBlob):
    """
    Check if two blobs overlaps or right next to each other.

    :param selfBlob: A :class:`pdb_eda.ccp4.DensityBlob` object.
    :param otherBlob: A :class:`pdb_eda.ccp4.DensityBlob` object.
    :return: :py:obj:`True` or :py:obj:`False`.
    """
    return True if any(-1 <= x[0] - y[0] <= 1 and -1 <= x[1] - y[1] <= 1 and -1 <= x[2] - y[2] <= 1 for x in selfBlob.crsList for y in otherBlob.crsList) else False


def sumOfAbs(array, float cutoff):
    return sum(abs(value) for value in array if abs(value) > cutoff)

import numpy as np
import scipy.spatial
dcutoff = np.sqrt(3)  ## the points are considered to be adjacent if one is in the one layer outer box with the other one in the center
def createCrsLists(crsList):
    """
    Calculates a list of crsLists from a given crsList.
    This is a preparation step for creating blobs.

    :param crsList: a crs list.
    :return: crsLists is a list of crsList.
    :rtype: A :py:obj:`list`
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
def createSymmetryAtoms(list atomList, rotationMats, orthoMat, list xs, list ys, list zs):
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
    """
    A wrapper class to the `BioPDB.atom` class, delegating all BioPDB atom class methods and data members except having its own symmetry and coordination.
    """

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

