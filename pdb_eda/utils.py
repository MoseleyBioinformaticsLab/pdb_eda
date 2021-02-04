# Compile using:
# $ python3 ccp4_utils_setup.py build_ext --inplace
# $ cp build/lib.linux-x86_64-3.6/pdb_eda/ccp4_utils.cpython-36m-x86_64-linux-gnu.so ccp4_utils.so

def testOverlap(selfBlob, otherBlob):
    """
    Check if two blobs overlaps or right next to each other.

    :param selfBlob: A :class:`pdb_eda.ccp4.DensityBlob` object.
    :param otherBlob: A :class:`pdb_eda.ccp4.DensityBlob` object.
    :return: :py:obj:`True` or :py:obj:`False`.
    """
    #if any(x in self.crsList for x in otherBlob.crsList):
    #    return True
    #if np.any(scipy.spatial.distance.cdist(np.matrix(self.crsList), np.matrix(otherBlob.crsList)) <= np.sqrt(3)):
    if any(-1 <= x[0] - y[0] <= 1 and -1 <= x[1] - y[1] <= 1 and -1 <= x[2] - y[2] <= 1 for x in selfBlob.crsList for y in otherBlob.crsList):
        return True
    else:
        return False


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
