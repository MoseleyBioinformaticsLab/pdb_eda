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
