"""
Cythonized Utilities (pdb_eda.cutils)
------------------------------------
"""

def testOverlap(selfBlob, otherBlob):
    """Check if two blobs overlaps or right next to each other.

    :param :class:`pdb_eda.ccp4.DensityBlob` selfBlob:
    :param :class:`pdb_eda.ccp4.DensityBlob` otherBlob:
    :return: bool
    :rtype: bool
    """
    return True if any(-1 <= x[0] - y[0] <= 1 and -1 <= x[1] - y[1] <= 1 and -1 <= x[2] - y[2] <= 1 for x in selfBlob.crsList for y in otherBlob.crsList) else False


def sumOfAbs(array, float cutoff):
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

    :param list crsList: a crs list.
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
    """Creates and returns a list of all symmetry atoms.

    :param list atomList:
    :param list rotationMats:
    :param list orthoMat:
    :param list xs:
    :param list ys:
    :param list zs:
    :return: allAtoms
    :rtype: list
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
        """`pdb_eda.densityAnalysis.symAtom` initializer.

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

    :param `pdb_eda.ccp4.DensityMatrix` densityMatrix:
    :param crsCoord: crs coordinates.
    :type: list or set
    :return: density
    :rtype: float
    """
    crsCoord = list(crsCoord)
    header = densityMatrix.header
    for ind in range(3):
        if crsCoord[ind] < 0 or crsCoord[ind] >= header.ncrs[ind]:
            crsCoord[ind] -= int(np.floor(crsCoord[ind] / header.crsInterval[ind]) * header.crsInterval[ind])

        if (header.ncrs[ind] <= crsCoord[ind] < header.crsInterval[ind]) or crsCoord[ind] < 0:
            return 0

    return densityMatrix.density[crsCoord[2], crsCoord[1], crsCoord[0]]

def testValidCrs(densityMatrix, crsCoord):
    """Get the density of a point.

    :param `pdb_eda.ccp4.DensityMatrix` densityMatrix:
    :param list crsCoord: crs coordinates.
    :return: bool
    :rtype: bool
    """
    crsCoord = list(crsCoord)
    header = densityMatrix.header
    for ind in range(3):
        if crsCoord[ind] < 0 or crsCoord[ind] >= header.ncrs[ind]:
            crsCoord[ind] -= int(np.floor(crsCoord[ind] / header.crsInterval[ind]) * header.crsInterval[ind])

        if (header.ncrs[ind] <= crsCoord[ind] < header.crsInterval[ind]) or crsCoord[ind] < 0:
            return False

    return True

def testValidCrsList(densityMatrix, crsList):
    """Get the density of a point.

    :param `pdb_eda.ccp4.DensityMatrix` densityMatrix:
    :param crsList: list of crs coordinates.
    :type: list or set
    :return: bool
    :rtype: bool
    """
    return not any(not testValidCrs(densityMatrix,crs) for crs in crsList)

def createFullCrsList(densityMatrix, float cutoff):
    """Returns full crs list for the density matrix.

    :param `pdb_eda.ccp4.DensityMatrix` densityMatrix:
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

cdef inline _testXyzWithinDistance(xyzCoord1, xyzCoord2, float distance):
    """Test whether two xyzCoords are within a certain distance.

    :param list xyzCoord1:
    :param list xyzCoord2:
    :param float radius:
    :return: bool
    :rtype: bool
    """
    return np.sqrt((xyzCoord2[0] - xyzCoord1[0])**2 + (xyzCoord2[1] - xyzCoord1[1])**2 + (xyzCoord2[2] - xyzCoord1[2])**2) <= distance

def getSphereCrsFromXyz(densityMatrix, xyzCoord, float radius, float densityCutoff=0):
    """Calculate a list of crs coordinates that within a given distance of a xyz point.

    :param `pdb_eda.ccp4.DensityMatrix` densityMatrix:
    :param list xyzCoord: xyz coordinates.
    :param float radius: the radius.
    :param float densityCutoff: a density cutoff for all the points wants to be included.
            Default 0 means include every point within the radius.
            If cutoff < 0, include only points with density < cutoff.
            If cutoff > 0, include only points with density > cutoff.

    :return: crsCoordList of crs coordinates
    :rtype: list
    """
    crsCoord = densityMatrix.header.xyz2crsCoord(xyzCoord)
    crsRadius = densityMatrix.header.xyz2crsCoord(densityMatrix.origin + [radius, radius, radius])
    crsCoordList = []
    for crs in itertools.product(range(crsCoord[0] - crsRadius[0]-1, crsCoord[0] + crsRadius[0]+1),
                                 range(crsCoord[1] - crsRadius[1]-1, crsCoord[1] + crsRadius[1]+1),
                                 range(crsCoord[2] - crsRadius[2]-1, crsCoord[2] + crsRadius[2]+1)):
        density = getPointDensityFromCrs(densityMatrix,crs)
        if ((0 < densityCutoff < density) or (density < densityCutoff < 0) or densityCutoff == 0) and _testXyzWithinDistance(xyzCoord,densityMatrix.header.crs2xyzCoord(crs),radius):
            crsCoordList.append(crs)

    return crsCoordList

def getSphereCrsFromXyzList(densityMatrix, xyzCoordList, float radius, float densityCutoff=0):
    """Calculate a list of crs coordinates that within a given distance from a list of xyz points.

    :param `pdb_eda.ccp4.DensityMatrix` densityMatrix:
    :param list xyzCoord: xyz coordinates.
    :param float radius: the radius.
    :param float densityCutoff: a density cutoff for all the points wants to be included.
            Default 0 means include every point within the radius.
            If cutoff < 0, include only points with density < cutoff.
            If cutoff > 0, include only points with density > cutoff.

    :return: crsCoordList of crs coordinates
    :rtype: set
    """
    return {tuple(crsCoord) for xyzCoord in xyzCoordList for crsCoord in getSphereCrsFromXyz(densityMatrix, xyzCoord, radius, densityCutoff)}

def testValidXyz(densityMatrix, xyzCoord, float radius):
    """Test whether all crs coordinates within a given distance of a xyzCoord is within the densityMatrix.

    :param `pdb_eda.ccp4.DensityMatrix` densityMatrix:
    :param list xyzCoord: xyz coordinates.
    :param float radius: the radius.
    :return: bool
    :rtype: bool
    """
    crsCoord = densityMatrix.header.xyz2crsCoord(xyzCoord)
    crsRadius = densityMatrix.header.xyz2crsCoord(densityMatrix.origin + [radius, radius, radius])
    return not any(not testValidCrs(densityMatrix, crs)
                   for crs in itertools.product(range(crsCoord[0] - crsRadius[0]-1, crsCoord[0] + crsRadius[0]+1),
                                                range(crsCoord[1] - crsRadius[1]-1, crsCoord[1] + crsRadius[1]+1),
                                                range(crsCoord[2] - crsRadius[2]-1, crsCoord[2] + crsRadius[2]+1))
                   if _testXyzWithinDistance(xyzCoord, densityMatrix.header.crs2xyzCoord(crs), radius))

def testValidXyzList(densityMatrix, xyzCoordList, float radius):
    """Test whether all crs coordinates within a given distance of a set of xyzCoords is within the densityMatrix.

    :param `pdb_eda.ccp4.DensityMatrix` densityMatrix:
    :param list xyzCoordList: list of xyz coordinates.
    :param float radius: the radius.
    :return: bool
    :rtype: bool
    """
    return not any(not testValidXyz(densityMatrix,xyzCoord,radius) for xyzCoord in xyzCoordList)