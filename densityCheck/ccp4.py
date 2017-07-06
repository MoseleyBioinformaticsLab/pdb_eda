# !/usr/bin/python3
"""
ccpp4.py
    Reads and parses the CCP4 format binary files and returns DensityMatrix objects.
    Format details of ccp4 can be found in http://www.ccp4.ac.uk/html/maplib.html

"""

import urllib.request
import struct
import numpy as np
import statistics
import warnings

urlPrefix = "http://www.ebi.ac.uk/pdbe/coordinates/files/"
urlSuffix = ".ccp4"


def readFromPDBID(pdbid, verbose=False):
    """RETURNS DensityMatrix object given the PARAMETER pdbid."""
    return readFromURL(urlPrefix + pdbid.lower() + urlSuffix, verbose)


def readFromURL(url, verbose=False):
    """RETURNS DensityMatrix object given the PARAMETER url."""
    with urllib.request.urlopen(url) as urlHandle:
        return parse(urlHandle, verbose)


def read(ccp4Filename, verbose=False):
    """RETURNS DensityMatrix object given the PARAMETER fileName."""
    with open(ccp4Filename, "rb") as fileHandle:
        return parse(fileHandle, verbose)


def parse(handle, verbose=False):
    """RETURNS DensityMatrix object given the PARAMETER file handle."""
    header = DensityHeader.fromFileHeader(handle.read(1024))
    endian = header.endian
    dataBuffer = handle.read()

    # Sanity check on file sizes
    if len(dataBuffer) != header.symmetryBytes + header.mapSize:
        assert header.symmetryBytes == 0 | len(
            dataBuffer) != header.mapSize, "Error: File contains suspicious symmetry records"
        assert header.mapSize == 0 | len(dataBuffer) != header.symmetryBytes, "Error: File contains no map data"
        assert len(dataBuffer) > header.symmetryBytes + header.mapSize, "Error: contains incomplete data"
        assert len(dataBuffer) < header.symmetryBytes + header.mapSize, "Error: File contains larger than expected data"

    assert header.xlength != 0.0 or header.ylength != 0.0 or header.zlength != 0.0, "Error: Cell dimensions are all 0, Map file will not align with other structures"

    if header.nintervalX == 0 & header.ncrs[0] > 0:
        header.nintervalX = header.ncrs[0] - 1
        if verbose: warnings.warn("Fixed number of X interval")
    if header.nintervalY == 0 & header.ncrs[1] > 0:
        header.nintervalY = header.ncrs[1] - 1
        if verbose: warnings.warn("Fixed number of Y interval")
    if header.nintervalZ == 0 & header.ncrs[2] > 0:
        header.nintervalZ = header.ncrs[2] - 1
        if verbose: warnings.warn("Fixed number of Z interval.")

    if header.col2xyz == 0 & header.row2xyz == 0 & header.sec2xyz == 0:
        header.col2xyz = 1
        header.row2xyz = 2
        header.sec2xyz = 3
        if verbose: warnings.warn("Mappings from column/row/section to xyz are all 0, set to 1, 2, 3 instead.")

    symmetry = dataBuffer[0:header.symmetryBytes]
    mapData = dataBuffer[header.symmetryBytes:len(dataBuffer)]

    numBytes = int(len(mapData) / 4)
    densities = struct.unpack(endian + numBytes * 'f', mapData)
    origin = header.origin

    # Calculate some statistics
    sigma = np.std(densities)
    mean = np.mean(densities)
    median = np.median(densities)
    mode = 0  # statistics.mode(densities)
    #print('mean, median, mode, sigma, header rmsd, difference of the last two: ', mean, median, mode, sigma, header.rmsd, sigma - header.rmsd)

    return DensityMatrix(header, origin, densities)


class DensityHeader(object):
    @classmethod
    def fromFileHeader(cls, fileHeader):
        """RETURNS DensityHeader object given the PARAMETER fileHeader."""

        # Test for endianness
        mode = int.from_bytes(fileHeader[12:16], byteorder='little')
        endian = '<' if 0 <= mode <= 6 else '>'

        # Header
        headerFormat = endian + 10 * 'i' + 6 * 'f' + 3 * 'i' + 3 * 'f' + 3 * 'i' + 27 * 'f' + 4 * 'c' + 'ifi'
        headerTuple = struct.unpack(headerFormat, fileHeader[:224])
        #print(headerTuple)
        labels = fileHeader[224:]  # Labels in header
        labels = labels.replace(b' ', b'')

        header = DensityHeader(headerTuple, labels, endian)
        return header

    def __init__(self, headerTuple, labels, endian):
        """
        Initialize the DensityHeader object, name data members accordingly and calculate some metric that will be used frequently
        PARAMS
            :param headerTuple: The ccp4 header information (excluding labels) in a tuple
            :param labels: The labels field in a ccp4 header.
            :param endian: the endianness of the file
        RETURNS
            DensityHeader object
        """
        self.ncrs = headerTuple[0:3]
        """
        Number of Columns    (fastest changing in map)
        Number of Rows
        Number of Sections   (slowest changing in map)
        """
        self.mode = headerTuple[3]
        self.endian = endian
        """
        Data type
            0 = envelope stored as signed bytes (from -128 lowest to 127 highest)
            1 = Image     stored as Integer*2
            2 = Image     stored as Reals
            3 = Transform stored as Complex Integer*2
            4 = Transform stored as Complex Reals
            5 == 0

            Note: Mode 2 is the normal mode used in the CCP4 programs. Other modes than 2 and 0
                may NOT WORK
        """
        self.crsStart = headerTuple[4:7]  # Number of first COLUMN, ROW, and SECTION in map
        self.nintervalX = headerTuple[7]  # Number of intervals along X
        self.nintervalY = headerTuple[8]  # Number of intervals along Y
        self.nintervalZ = headerTuple[9]  # Number of intervals along Z
        self.xlength = headerTuple[10]  # Cell Dimensions (Angstroms)
        self.ylength = headerTuple[11]  # ''
        self.zlength = headerTuple[12]  # ''
        self.alpha = headerTuple[13]  # Cell Angles     (Degrees)
        self.beta = headerTuple[14]  # ''
        self.gamma = headerTuple[15]  # ''
        self.col2xyz = headerTuple[16]  # Which axis corresponds to Cols.  (1,2,3 for X,Y,Z)
        self.row2xyz = headerTuple[17]  # Which axis corresponds to Rows   (1,2,3 for X,Y,Z)
        self.sec2xyz = headerTuple[18]  # Which axis corresponds to Sects. (1,2,3 for X,Y,Z)
        self.densityMin = headerTuple[19]  # Minimum density value
        self.densityMax = headerTuple[20]  # Maximum density value
        self.densityMean = headerTuple[21]  # Mean    density value    (Average)
        self.spaceGroup = headerTuple[22]  # Space group number
        self.symmetryBytes = headerTuple[23]  # Number of bytes used for storing symmetry operators
        self.skewFlag = headerTuple[24]  # Flag for skew transformation, =0 none, =1 if foll
        self.skewMat = headerTuple[25:34]  # Skew matrix S (in order S11, S12, S13, S21 etc) if LSKFLG .ne. 0.
        self.skewTrans = headerTuple[34:37]
        """
        Skew translation t if LSKFLG .ne. 0.
                    Skew transformation is from standard orthogonal
                    coordinate frame (as used for atoms) to orthogonal
                    map frame, as: Xo(map) = S * (Xo(atoms) - t)
        """
        self.futureUse = headerTuple[37:49]
        """
        (some of these are used by the MSUBSX routines in MAPBRICK, MAPCONT and
        FRODO) (all set to zero by default)
        """
        self.originEM = headerTuple[49:52]
        """
        Use ORIGIN records rather than old crsStart records as in http://www2.mrc-lmb.cam.ac.uk/image2000.html
        The ORIGIN field is only used by the EM community, and has undefined meaning for non-orthogonal maps and/or
        non-cubic voxels, etc.
        """
        self.mapChar = headerTuple[52:56]  # Character string 'MAP ' to identify file type
        self.machineStamp = headerTuple[56]  # Machine stamp indicating the machine type which wrote file
        self.rmsd = headerTuple[57]  # Rms deviation of map from mean density
        self.nLabel = headerTuple[58]  # Number of labels being used
        self.labels = labels

        self.mapSize = self.ncrs[0] * self.ncrs[1] * self.ncrs[2] * 4
        self.xyzLength = [self.xlength, self.ylength, self.zlength]
        self.xyzInterval = [self.nintervalX, self.nintervalY, self.nintervalZ]
        self.gridLength = [x/y for x, y in zip(self.xyzLength, self.xyzInterval)]

        indices = [0, 0, 0]
        indices[self.col2xyz - 1] = 0
        indices[self.row2xyz - 1] = 1
        indices[self.sec2xyz - 1] = 2
        self.map2xyz = indices
        self.map2crs = [self.col2xyz - 1, self.row2xyz - 1, self.sec2xyz - 1]

        alpha = np.pi / 180 * self.alpha
        beta = np.pi / 180 * self.beta
        gamma = np.pi / 180 * self.gamma
        self.unitVolume = self.xlength * self.ylength * self.zlength / self.nintervalX / self.nintervalY / self.nintervalZ * \
                          np.sqrt(1 - np.cos(alpha) ** 2 - np.cos(beta) ** 2 - np.cos(gamma) ** 2 + 2 * np.cos(alpha) * np.cos(beta) * np.cos(gamma))

        ## A resuable part in the cell volumn calculation
        temp = np.sqrt(1 - np.cos(alpha) ** 2 - np.cos(beta) ** 2 - np.cos(gamma) ** 2 + 2 * np.cos(alpha) * np.cos(beta) * np.cos(gamma))

        self.orthoMat = [[self.xlength, self.ylength * np.cos(gamma), self.zlength * np.cos(beta)],
                         [0, self.ylength * np.sin(gamma), self.zlength * (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma)],
                         [0, 0, self.zlength * temp / np.sin(gamma)]]

        self.deOrthoMat = np.linalg.inv(self.orthoMat)
        self.deOrthoMat[abs(self.deOrthoMat) < 1e-10] = 0.0

        """
        self.deOrthoMat = [[1/self.xlength, - np.cos(gamma) / np.sin(gamma) / self.xlength,
                            (np.cos(gamma) * (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma) - np.cos(beta) * np.sin(gamma)) / self.xlength / temp],
                           [0, 1/np.sin(gamma)/self.ylength, - (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma) / self.ylength / temp],
                           [0, 0, np.sin(gamma) / self.zlength / temp]]
        """
        self.origin = self._calculateOrigin()


    def _calculateOrigin(self):
        """
        Calculate and RETURNS the xyz coordinates from the header information with no extra PARAMETER
        """
        # Orthogonalization matrix for calculation between fractional coordinates and orthogonal coordinates
        # Formula based on 'Biomolecular Crystallography' by Bernhard Rupp, p233

        if self.futureUse[-3] == 0.0 and self.futureUse[-2] == 0.0 and self.futureUse[-1] == 0.0:
            origin = np.dot(self.orthoMat, [self.crsStart[self.map2xyz[i]] / self.xyzInterval[i] for i in range(3)])
        else:
            origin = [self.originEM[i] for i in range(3)]

        return origin


    def xyz2crsCoord(self, xyzCoord):
        """
        Convert the xyz coordinates into crs coordinates
        PARAMS
            xyzCoord: xyz coordinates
        RETURNS
            crs coordinates
        """
        if self.alpha == self.beta == self.gamma == 90:
            crsGridPos = [int(round((xyzCoord[i] - self.origin[i]) / self.gridLength[i])) for i in range(3)]
        else:
            fraction = np.dot(self.deOrthoMat, xyzCoord)
            crsGridPos = [int(fraction[i] * self.xyzInterval[i]) - self.crsStart[self.map2xyz[i]] for i in range(3)]
        return [crsGridPos[self.map2crs[i]] for i in range(3)]

    def crs2xyzCoord(self, crsCoord):
        """
        Convert the crs coordinates into xyz coordinates
        PARAMS
            crsCoord: crs coordinates
        RETURNS
            xyz coordinates
        """
        if self.alpha == self.beta == self.gamma == 90:
            return [int(crsCoord[self.map2xyz[i]]) * self.gridLength[i] + self.origin[i] for i in range(3)]
        else:
            return np.dot(self.orthoMat, [(crsCoord[self.map2xyz[i]] + self.crsStart[self.map2xyz[i]]) / self.xyzInterval[i] for i in range(3)])


class DensityMatrix:
    def __init__(self, header, origin, density):
        """
        Initialize a DensityMatrix object
        PARAMS
            :param header: the DensityHeader object of the density matrix
            :param origin: the xyz coordinates of the origin of the first number of the density data
            :param density: the density data as a 1-d list
        RETURNs
            DensityMatrix object
        """
        self.header = header
        self.origin = origin
        self.densityArray = density
        self.density = np.array(density).reshape(header.ncrs[2], header.ncrs[1], header.ncrs[0])


    def validCRS(self, crsCoord):
        """Check if the crs coordinate is valid (within the range of data) given PARAMETER crs coordinate"""
        for ind in range(3):
            crsInterval = self.header.xyzInterval[self.header.map2crs[ind]]

            if crsCoord[ind] < 0 or crsCoord[ind] >= self.header.ncrs[ind]:
                n = np.floor(crsCoord[ind] / crsInterval)
                crsCoord[ind] -= int(n * crsInterval)

            if self.header.ncrs[ind] <= crsCoord[ind] < crsInterval:
                return False

        return True

    def getPointDensityFromCrs(self, crsCoord):
        """RETURNS the density of a point given PARAMETER crs coordinate"""
        if not self.validCRS(crsCoord):
            # message = "No data available in " + str(ind + 1) + " axis"
            # warnings.warn(message)
            return 0

        return self.density[crsCoord[2], crsCoord[1], crsCoord[0]]

    def getPointDensityFromXyz(self, xyzCoord):
        """RETURNS the density of a point given PARAMETER xyz coordinate"""
        crsCoord = self.header.xyz2crsCoord(xyzCoord)
        # print("crs grid: ", crsCoord)

        return self.getPointDensityFromCrs(crsCoord)

    def getSphereCrsFromXyz(self, xyzCoord, radius, densityCutoff=0):
        """
        RETURNS
            A list of crs coordinates that meet all criteria given
        PARAMS
            xyzCoord: xyz coordinates
            radius: the radius
            densityCutoff: a density cutoff for all the points wants to be included.
                    Default 0 means include every point within the radius.
                    If cutoff < 0, include only points with density < cutoff.
                    If cutoff > 0, include only points with density > cutoff.
        """
        crsCoord = self.header.xyz2crsCoord(xyzCoord)

        xyzRadius = [np.ceil(radius / self.header.gridLength[i]) for i in range(3)]
        crsRadius = [int(x) for x in [xyzRadius[self.header.map2crs[y]] for y in range(3)]]

        # print('grid positions', crsCoord)
        # print("crs radius:", crsRadius)
        # print('cutoff: ', densityCutoff)
        crsCoordList = []
        for cInd in range(-crsRadius[0], crsRadius[0]+1):
            for rInd in range(-crsRadius[1], crsRadius[1]+1):
                for sInd in range(- crsRadius[2], crsRadius[2]+1):
                    crs = [x + y for x, y in zip(crsCoord, [cInd, rInd, sInd])]
                    if 0 < densityCutoff < self.getPointDensityFromCrs(crs) or self.getPointDensityFromCrs(crs) < densityCutoff < 0 or densityCutoff == 0:
                        if self.header.alpha == self.header.beta == self.header.gamma == 90:
                            if cInd**2 / crsRadius[0]**2 + rInd**2 / crsRadius[1]**2 + sInd**2 / crsRadius[2]**2 <= 1:
                                crsCoordList.append(crs)
                        else:
                            xyz = self.header.crs2xyzCoord(crs)
                            dist = np.sqrt((xyz[0] - xyzCoord[0])**2 + (xyz[1] - xyzCoord[1])**2 + (xyz[2] - xyzCoord[2])**2)
                            if dist < radius:
                                crsCoordList.append(crs)

        # print('crs grids: ', crsCoordList)

        return crsCoordList

    def getTotalDensityFromXyz(self, xyzCoord, radius, densityCutoff=0):
        """
        RETURNS a total density given
        PARAMS
            xyzCoord: xyz coordinates
            radius: the radius
            densityCutoff: a density cutoff for all the points wants to be included.
                    Default 0 means include every point within the radius.
                    If cutoff < 0, include only points with density < cutoff.
                    If cutoff > 0, include only points with density > cutoff.
        """
        crsCoordList = self.getSphereCrsFromXyz(xyzCoord, radius, densityCutoff)

        totalDensity = 0
        for crs in crsCoordList:
            totalDensity += self.getPointDensityFromCrs(crs)

        return totalDensity

    def findAberrantBlobs(self, xyzCoord, radius, densityCutoff=0):
        """
        Within given range, find and aggregates all neighbouring aberrant points into blobs (red/green meshes) and
        RETURNS
            a list of aberrant blobs described by their xyz centroid, total density, and volume.
        PARAMS
             xyzCoord: the xyz coordinates of a center
             radius: the radius
             densityCutoff: a density cutoff for all the points wants to be included.
                    Default 0 means include every point within the radius.
                    If cutoff < 0, include only points with density < cutoff.
                    If cutoff > 0, include only points with density > cutoff.
        """
        crsCoordList = self.getSphereCrsFromXyz(xyzCoord, radius, densityCutoff)

        adjSetList = []
        usedIndex = set()
        for ind in range(len(crsCoordList)):
            if ind in usedIndex:
                continue

            oldSet = set()
            newSet = {ind}
            while True:
                diffSet = newSet - oldSet
                adjacentSet = {x for x in range(len(crsCoordList)) for y in diffSet if -1 <= (crsCoordList[y][0] - crsCoordList[x][0]) <= 1 and -1 <= (crsCoordList[y][1] - crsCoordList[x][1]) <= 1 and -1 <= (crsCoordList[y][2] - crsCoordList[x][2]) <= 1}

                oldSet = newSet.copy()
                newSet.update(adjacentSet)
                if oldSet == newSet:
                    break

            usedIndex.update(newSet)
            adjSetList.append(newSet)

        # print("New function: ", adjSetList)

        blobs = []
        for adjacentSet in adjSetList:
            coords = [crsCoordList[x] for x in adjacentSet]

            blob = DensityBlob.fromCrsList(coords, self.header, self.density)
            blobs.append(blob)

        return blobs


class DensityBlob:
    def __init__(self, centroid, totalDensity, volume, crsList, header):
        """
        Initialize a DensityBlob object
        PARAMS
            :param centroid: the centroid of the blob
            :param totalDensity: the totalDensity of the blob
            :param volume: the volume of the blob, number of density units * unit volumes
            :param crsList: the crs list of the blob
            :param header: the header of the ccp4 file
        RETURNs
            DensityBlob object
        """
        self.centroid = centroid
        self.totalDensity = totalDensity
        self.volume = volume
        self.crsList = crsList
        self.header = header


    @staticmethod
    def fromCrsList(crsList, header, densityMat):
        """
        The creator of a DensityBlob object
        PARAMS
            crsList: the crs list of the blob
            header: the header of the ccp4 file
            densityMat: the 3-d density matrix that just passed in for calculating centroid etc, so the object does not have to have a density list data member
        RETURNS
            DensityBlob object
        """
        weights = [0, 0, 0]
        totalDen = 0
        for i, point in enumerate(crsList):
            density = densityMat[point[2], point[1], point[0]]
            pointXYZ = header.crs2xyzCoord(point)
            weights = [weights[i] + density * pointXYZ[i] for i in range(3)]
            totalDen += density

        centroidXYZ = [weight / totalDen for weight in weights]
        return DensityBlob(centroidXYZ, totalDen, header.unitVolume * len(crsList), crsList, header)


    def __eq__(self, other):
        """
        OVERRIDE the '==' operator of DensityBlob object
        """
        if abs(self.volume - other.volume) >= 1e-6: return False
        if abs(self.totalDensity - other.totalDensity) >= 1e-6: return False
        for i in range(0, 3):
            if abs(self.centroid[i] - other.centroid[i]) >= 1e-6: return False

        return True

    def testOverlap(self, otherBlob):
        """
        RETURN true if two blobs overlaps or right next to each other
        """
        if any(x in self.crsList for x in otherBlob.crsList):
            return True
        elif any(-1 <= x[0] - y[0] <= 1 and -1 <= x[1] - y[1] <= 1 and -1 <= x[2] - y[2] <= 1 for x in self.crsList for y in otherBlob.crsList):
            return True
        else:
            return False

    def merge(self, otherBlob, densityMat):
        """
        UPDATE the original blob given PARAM another blob
        """
        combinedList = self.crsList + [x for x in otherBlob.crsList if x not in self.crsList]
        newBlob = DensityBlob.fromCrsList(combinedList, self.header, densityMat)

        self.centroid = newBlob.centroid
        self.totalDensity = newBlob.totalDensity
        self.volume = newBlob.volume
        self.crsList = newBlob.crsList
        self.header = newBlob.header



