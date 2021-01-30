def testOverlap(selfBlob, otherBlob):
    """
    Check if two blobs overlaps or right next to each other.

    :return: :py:obj:`True` or :py:obj:`False`.
    """
    #if any(x in self.crsList for x in otherBlob.crsList):
    #    return True
    #if np.any(scipy.spatial.distance.cdist(np.matrix(self.crsList), np.matrix(otherBlob.crsList)) <= np.sqrt(3)):
    if any(-1 <= x[0] - y[0] <= 1 and -1 <= x[1] - y[1] <= 1 and -1 <= x[2] - y[2] <= 1 for x in selfBlob.crsList for y in otherBlob.crsList):
        return True
    else:
        return False

