# !/usr/bin/python3
"""
CCP4 Parser (pdb_eda.ccp4)
-------------------------------------------------------

This module provides methods to read and parse the CCP4 format files, returning ccp4 objects.
Format details of ccp4 can be found in http://www.ccp4.ac.uk/html/maplib.html.
"""

import warnings
import struct

import urllib.request
import numpy as np

try:
    from . import cutils as utils
except ImportError:
    from . import utils

urlPrefix = "http://www.ebi.ac.uk/pdbe/coordinates/files/"
urlSuffix = ".ccp4"


def readFromPDBID(pdbid, verbose=False):
    """Creates :class:`pdb_eda.ccp4.DensityMatrix` object.

    :param pdbid: PDB entry ID.
    :type pdbid: :py:class:`str`
    :param verbose: verbose mode, defaults to :py:obj:`False`
    :type verbose: :py:class:`bool`

    :return: densityMatrix
    :rtype: :class:`pdb_eda.ccp4.DensityMatrix`
    """
    return readFromURL(urlPrefix + pdbid.lower() + urlSuffix, pdbid, verbose)


def readFromURL(url, pdbid=None, verbose=False):
    """Creates :class:`pdb_eda.ccp4.DensityMatrix` object.

    :param url:
    :type url: :py:class:`str`
    :param pdbid: PDB entry ID.
    :type pdbid: :py:class:`str`, optional
    :param verbose: verbose mode, defaults to :py:obj:`False`
    :type verbose: :py:class:`bool`

    :return: densityMatrix
    :rtype: :class:`pdb_eda.ccp4.DensityMatrix`
    """
    if not pdbid:
        pdbid = url
    with urllib.request.urlopen(url) as urlHandle:
        return parse(urlHandle, pdbid, verbose)


def read(ccp4Filename, pdbid=None, verbose=False):
    """Creates :class:`pdb_eda.ccp4.DensityMatrix` object.

    :param ccp4Filename: .ccp4 filename including path.
    :type ccp4Filename: :py:class:`str`
    :param pdbid: PDB entry ID.
    :type pdbid: :py:class:`str`, optional
    :param verbose: verbose mode, defaults to :py:obj:`False`
    :type verbose: :py:class:`bool`

    :return: densityMatrix
    :rtype: :class:`pdb_eda.ccp4.DensityMatrix`
    """
    if not pdbid:
        pdbid = ccp4Filename
    with open(ccp4Filename, "rb") as fileHandle:
        return parse(fileHandle, pdbid, verbose)


def parse(handle, pdbid, verbose=False):
    """Creates :class:`pdb_eda.ccp4.DensityMatrix` object.

    :param handle: an I/O handle for .ccp4 file.
    :type handle: :class:`io.IOBase`
    :param pdbid: PDB entry ID.
    :type pdbid: :py:class:`str`
    :param verbose: verbose mode, defaults to :py:obj:`False`
    :type verbose: :py:class:`bool`

    :return: densityMatrix
    :rtype: :class:`pdb_eda.ccp4.DensityMatrix`
    """
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

    header.symmetry = dataBuffer[0:header.symmetryBytes]
    mapData = dataBuffer[header.symmetryBytes:len(dataBuffer)]

    numBytes = int(len(mapData) / 4)
    densities = struct.unpack(endian + numBytes * 'f', mapData)
    origin = header.origin

    return DensityMatrix(header, origin, densities, pdbid)


class DensityHeader(object):
    """:class:`pdb_eda.ccp4.DensityHeader` that stores information about ccp4 header."""

    @classmethod
    def fromFileHeader(cls, fileHeader):
        """RETURNS :class:`pdb_eda.ccp4.DensityHeader` object given the fileHeader.

        :param fileHeader: ccp4 file header.
        :type fileHeader: :py:class:`bytes`

        :return: densityHeader
        :rtype: :class:`pdb_eda.ccp4.DensityHeader`
        """

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
        """Initialize the DensityHeader object, assign values to data members accordingly, and calculate some metrics that will be used frequently.

        :param headerTuple: The ccp4 header information (excluding labels) in a tuple.
        :type headerTuple: :py:class:`tuple`
        :param labels: The labels field in a ccp4 header.
        :type labels: :py:class:`bytes`
        :param endian: The endianness of the file.
        :type endian: :py:class:`str`
        """
        self.ncrs = headerTuple[0:3]
        #Number of Columns    (fastest changing in map)
        #Number of Rows
        #Number of Sections   (slowest changing in map)

        self.mode = headerTuple[3]
        self.endian = endian
        #Data type
        #    0 = envelope stored as signed bytes (from -128 lowest to 127 highest)
        #    1 = Image     stored as Integer*2
        #    2 = Image     stored as Reals
        #    3 = Transform stored as Complex Integer*2
        #    4 = Transform stored as Complex Reals
        #    5 == 0

        #    Note: Mode 2 is the normal mode used in the CCP4 programs. Other modes than 2 and 0
        #        may NOT WORK

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
        #Skew translation t if LSKFLG .ne. 0.
        #            Skew transformation is from standard orthogonal
        #            coordinate frame (as used for atoms) to orthogonal
        #            map frame, as: Xo(map) = S * (Xo(atoms) - t)

        self.futureUse = headerTuple[37:49]
        #(some of these are used by the MSUBSX routines in MAPBRICK, MAPCONT and FRODO) (all set to zero by default)
        self.originEM = headerTuple[49:52]
        #Use ORIGIN records rather than old crsStart records as in http://www2.mrc-lmb.cam.ac.uk/image2000.html
        #The ORIGIN field is only used by the EM community, and has undefined meaning for non-orthogonal maps and/or
        #non-cubic voxels, etc.

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

        self.crsInterval = [self.xyzInterval[self.map2crs[ind]] for ind in range(3)]


        alpha = np.pi / 180 * self.alpha
        beta = np.pi / 180 * self.beta
        gamma = np.pi / 180 * self.gamma
        self.unitVolume = self.xlength * self.ylength * self.zlength / self.nintervalX / self.nintervalY / self.nintervalZ * \
                          np.sqrt(1 - np.cos(alpha) ** 2 - np.cos(beta) ** 2 - np.cos(gamma) ** 2 + 2 * np.cos(alpha) * np.cos(beta) * np.cos(gamma))

        ## A reusable part in the cell volumn calculation
        temp = np.sqrt(1 - np.cos(alpha) ** 2 - np.cos(beta) ** 2 - np.cos(gamma) ** 2 + 2 * np.cos(alpha) * np.cos(beta) * np.cos(gamma))
        self.orthoMat = [[self.xlength, self.ylength * np.cos(gamma), self.zlength * np.cos(beta)],
                         [0, self.ylength * np.sin(gamma), self.zlength * (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma)],
                         [0, 0, self.zlength * temp / np.sin(gamma)]]

        self.deOrthoMat = np.linalg.inv(self.orthoMat)
        self.deOrthoMat[abs(self.deOrthoMat) < 1e-10] = 0.0

        #self.deOrthoMat = [[1/self.xlength, - np.cos(gamma) / np.sin(gamma) / self.xlength,
        #                    (np.cos(gamma) * (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma) - np.cos(beta) * np.sin(gamma)) / self.xlength / temp],
        #                   [0, 1/np.sin(gamma)/self.ylength, - (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma) / self.ylength / temp],
        #                   [0, 0, np.sin(gamma) / self.zlength / temp]]

        self.origin = self._calculateOrigin()

        ncrs = [i for i in self.ncrs]
        if self.xyzInterval[self.col2xyz - 1] < self.ncrs[0]:
            ncrs[0] = self.xyzInterval[self.col2xyz - 1]
        if self.xyzInterval[self.row2xyz - 1] < self.ncrs[1]:
            ncrs[1] = self.xyzInterval[self.row2xyz - 1]
        if self.xyzInterval[self.sec2xyz - 1] < self.ncrs[2]:
            ncrs[2] = self.xyzInterval[self.sec2xyz - 1]
        self.uniqueNcrs = ncrs


    def _calculateOrigin(self):
        """Calculate the xyz coordinates from the header information.

        :return: xyz coordinates.
        :rtype: :py:class:`list` of :py:class:`float`.
        """
        # Orthogonalization matrix for calculation between fractional coordinates and orthogonal coordinates
        # Formula based on 'Biomolecular Crystallography' by Bernhard Rupp, p233

        if self.futureUse[-3] == 0.0 and self.futureUse[-2] == 0.0 and self.futureUse[-1] == 0.0:
            origin = np.dot(self.orthoMat, [self.crsStart[self.map2xyz[i]] / self.xyzInterval[i] for i in range(3)])
        else:
            origin = [self.originEM[i] for i in range(3)]

        return origin

    def xyz2crsCoord(self, xyzCoord):
        """Convert the xyz coordinates into crs coordinates.

        :param xyzCoord: xyz coordinates.
        :type xyzCoord: :py:class:`list` of :py:class:`float`

        :return: crs coordinates.
        :rtype: A :py:class:`list` of :py:class:`int`.
        """
        if self.alpha == self.beta == self.gamma == 90:
            crsGridPos = [int(round((xyzCoord[i] - self.origin[i]) / self.gridLength[i])) for i in range(3)]
        else:
            fraction = np.dot(self.deOrthoMat, xyzCoord)
            crsGridPos = [int(round(fraction[i] * self.xyzInterval[i])) - self.crsStart[self.map2xyz[i]] for i in range(3)]
        return [crsGridPos[self.map2crs[i]] for i in range(3)]

    def crs2xyzCoord(self, crsCoord):
        """Convert the crs coordinates into xyz coordinates.

        :param crsCoord: crs coordinates.
        :type crsCoord: A :py:obj:`list` of :py:class:`int`

        :return: xyz coordinates.
        :rtype: A :py:class:`list` of :py:class:`float`.
        """
        if self.alpha == self.beta == self.gamma == 90:
            return [crsCoord[self.map2xyz[i]] * self.gridLength[i] + self.origin[i] for i in range(3)]
        else:
            return np.dot(self.orthoMat, [(crsCoord[self.map2xyz[i]] + self.crsStart[self.map2xyz[i]]) / self.xyzInterval[i] for i in range(3)])


class DensityMatrix:
    """:class:`pdb_eda.ccp4.DensityMatrix` that stores data and methods of a ccp4 file."""

    def __init__(self, header, origin, density, pdbid):
        """Initialize the :class:`pdb_eda.ccp4.DensityMatrix` object.

        :param header:
        :type header: :class:`pdb_eda.ccp4.DensityHeader`
        :param origin: the xyz coordinates of the origin of the first number of the density data.
        :type origin: :py:class:`list`
        :param density: the density data as a 1-d list.
        :type density: :py:class:`tuple`
        :param pdbid: PDB entry ID
        :type pdbid: :py:class:`str`
        """
        self.pdbid = pdbid
        self.header = header
        self.origin = origin
        self.densityArray = density
        self.density = np.array(density).reshape(header.ncrs[2], header.ncrs[1], header.ncrs[0])
        self._meanDensity = None
        self._stdDensity = None
        self._totalAbsDensity = {}

    @property
    def meanDensity(self):
        """Returns mean of the density.

        :return: mean
        :rtype: :py:class:`float`
        """
        if self._meanDensity == None:
            self._meanDensity = np.mean(self.densityArray)
        return self._meanDensity

    @property
    def stdDensity(self):
        """Returns the standard deviation of the density.

        :return: std
        :rtype: :py:class:`float`
        """
        if self._stdDensity == None:
            self._stdDensity = np.std(self.densityArray)
        return self._stdDensity

    def getTotalAbsDensity(self, densityCutoff):
        """Returns total absolute Density above a densityCutoff

        :param densityCutoff:
        :type densityCutoff: :py:class:`float`

        :return: totalAbsDensity
        :rtype: :py:class:`float`
        """
        if densityCutoff not in self._totalAbsDensity:
            self._totalAbsDensity[densityCutoff] = utils.sumOfAbs(self.densityArray, densityCutoff)
        return self._totalAbsDensity[densityCutoff]

    def getPointDensityFromCrs(self, crsCoord):
        """Get the density of a point.

        :param crsCoord: crs coordinates.
        :type crsCoord: A :py:class:`list` of :py:class:`int`

        :return: pointDensity
        :rtype: :py:class:`float`
        """
        return utils.getPointDensityFromCrs(self,crsCoord)

    def getPointDensityFromXyz(self, xyzCoord):
        """Get the density of a point.

        :param xyzCoord: xyz coordinates.
        :type xyzCoord: A :py:class:`list` of :py:class:`float`

        :return: pointDensity
        :rtype: :py:class:`float`
        """
        return utils.getPointDensityFromCrs(self, self.header.xyz2crsCoord(xyzCoord))

    def getSphereCrsFromXyz(self, xyzCoord, radius, densityCutoff=0):
        """Calculate a list of crs coordinates that within a given distance of a point.

        :param xyzCoord: xyz coordinates.
        :type xyzCoord: A :py:class:`list` of :py:class:`float`
        :param radius:
        :type radius: :py:class:`float`
        :param densityCutoff: a density cutoff for all the points wants to be included, defaults to 0
                Default 0 means include every point within the radius.
                If cutoff < 0, include only points with density < cutoff.
                If cutoff > 0, include only points with density > cutoff.
        :type densityCutoff: :py:class:`float`

        :return: crsList
        :rtype: :py:class:`list`
        """
        return utils.getSphereCrsFromXyz(self,xyzCoord,radius,densityCutoff)

    def getTotalDensityFromXyz(self, xyzCoord, radius, densityCutoff=0):
        """Calculate the total density of a sphere.

        :param xyzCoord: xyz coordinates.
        :type xyzCoord: :py:class:`list` of :py:class:`float`
        :param radius:
        :type radius: :py:class:`float`
        :param densityCutoff: a density cutoff for all the points to include, defaults to 0
                Default 0 means include every point within the radius.
                If cutoff < 0, include only points with density < cutoff.
                If cutoff > 0, include only points with density > cutoff.
        :type densityCutoff: :py:class:`float`

        :return: density
        :rtype: :py:class:`float`
        """
        crsCoordList = utils.getSphereCrsFromXyz(self, xyzCoord, radius, densityCutoff)
        return sum(utils.getPointDensityFromCrs(self, crs) for crs in crsCoordList)

    def findAberrantBlobs(self, xyzCoords, radius, densityCutoff=0):
        """Within a given radius, find and aggregate all neighbouring aberrant points into blobs (red/green meshes).

        :param xyzCoords: single xyz coordinate or a list of xyz coordinates.
        :type xyzCoords: :py:class:`list`
        :param radius:
        :type radius: :py:class:`float`
        :param densityCutoff: A density cutoff for all the points wants to be included, defaults to 0
                Default 0 means include every point within the radius.
                If cutoff < 0, include only points with density < cutoff.
                If cutoff > 0, include only points with density > cutoff.
        :type densityCutoff: :py:class:`float`

        :return: blobList  of aberrant blobs described by their xyz centroid, total density, and volume.
        :rtype: :py:class:`list` of :class:`pdb_eda.ccp4.DensityBlob` objects.
        """
        if not isinstance(xyzCoords[0], (np.floating, float)): # test if xyzCoords is a single xyzCoord or a list of them.
            if len(xyzCoords) > 1:
                crsCoordList = list(utils.getSphereCrsFromXyzList(self, xyzCoords, radius, densityCutoff))
            else:
                crsCoordList = utils.getSphereCrsFromXyz(self, xyzCoords[0], radius, densityCutoff)
        else:
            crsCoordList = utils.getSphereCrsFromXyz(self, xyzCoords, radius, densityCutoff)

        return self.createBlobList(crsCoordList)

    def createFullBlobList(self, cutoff):
        """Aggregate the density map into positive (green or blue) or negative (red) blobs.

        :param cutoff: density cutoff to use to filter voxels.
        :type cutoff: :py:class:`float`

        :return blobList: list of DensityBlobs
        :rtype: :py:class:`list` of :class:`pdb_eda.ccp4.DensityBlob` objects.
        """
        crsList = utils.createFullCrsList(self, cutoff)
        return self.createBlobList(crsList) if crsList != None else None

    def createBlobList(self, crsList):
        """Calculates a list of blobs from a given crsList.

        :param crsList: a crs list.
        :type crsList: :py:class:`list`, :py:class:`set`

        :return: blobList
        :rtype: :py:class:`list` of :class:`pdb_eda.ccp4.DensityBlob` objects.
        """
        crsLists = utils.createCrsLists(crsList)
        return [ DensityBlob.fromCrsList(crs_list, self) for crs_list in crsLists ]


class DensityBlob:
    """:class:`pdb_eda.ccp4.DensityBlob` that stores data and methods of a electron density blob."""

    def __init__(self, centroid, coordCenter, totalDensity, volume, crsList, densityMatrix, atoms=None):
        """Initialize a :class:`pdb_eda.ccp4.DensityBlob` object.

        :param centroid: the centroid of the blob.
        :type centroid: :py:class:`list`
        :param totalDensity: the totalDensity of the blob.
        :type totalDensity: :py:class:`float`
        :param volume: the volume of the blob = number of density units * unit volumes.
        :type volume: :py:class:`float`
        :param crsList: the crs list of the blob.
        :type crsList: :py:class:`list`
        :param densityMatrix: the entire density map that the blob belongs to.
        :type densityMatrix: `pdb_eda.ccp4.DensityMatrix`
        :param atoms: list of atoms for the blob.
        :type atoms: :py:class:`list`, optional

        :return: densityBlob
        :rtype: :class:`pdb_eda.ccp4.DensityBlob`
        """
        self.centroid = centroid
        self.coordCenter = coordCenter
        self.totalDensity = totalDensity
        self.volume = volume
        self.crsList = {tuple(crs) for crs in crsList}
        self.densityMatrix = densityMatrix
        self.atoms = [] if not atoms else atoms

    @property
    def validCrs(self):
        return utils.testValidCrsList(self.densityMatrix, self.crsList)

    @staticmethod
    def fromCrsList(crsList, densityMatrix):
        """The creator of a A :class:`pdb_eda.ccp4.DensityBlob` object.

        :param crsList: the crs list of the blob.
        :type crsList: :py:class:`list`
        :param densityMatrix: the 3-d density matrix to use for calculating centroid etc, so the object does not have to have a density list data member.
        :type densityMatrix: :class:`pdb_eda.ccp4.DensityMatrix`

        :return: densityBlob
        :rtype: :class:`pdb_eda.ccp4.DensityBlob`
        """
        weights = [0, 0, 0]
        totalDensity = 0
        for i, point in enumerate(crsList):
            density = utils.getPointDensityFromCrs(densityMatrix, point)
            pointXYZ = densityMatrix.header.crs2xyzCoord(point)
            weights = [weights[i] + density * pointXYZ[i] for i in range(3)]
            totalDensity += density

        centroidXYZ = [weight / totalDensity for weight in weights]
        npoints = len(crsList)
        coordCenter = [sum(k) / npoints for k in zip(*[densityMatrix.header.crs2xyzCoord(crs) for crs in crsList])]
        return DensityBlob(centroidXYZ, coordCenter, totalDensity, densityMatrix.header.unitVolume * len(crsList), crsList, densityMatrix)


    def __eq__(self, otherBlob):
        """Check if two blobs are the same, and overwrite the '==' operator for the :class:`pdb_eda.ccp4.DensityBlob` object.

        :param otherBlob:
        :type otherBlob: :class:`pdb_eda.ccp4.DensityBlob`

        :return: bool
        :rtype: :py:class`bool`
        """
        if abs(self.volume - otherBlob.volume) >= 1e-6: return False
        if abs(self.totalDensity - otherBlob.totalDensity) >= 1e-6: return False
        for i in range(0, 3):
            if abs(self.centroid[i] - otherBlob.centroid[i]) >= 1e-6: return False

        return True

    def testOverlap(self, otherBlob):
        """Check if two blobs overlaps or right next to each other.

        :param otherBlob:
        :type otherBlob: :class:`pdb_eda.ccp4.DensityBlob`

        :return: bool
        :rtype: :py:class`bool`
        """
        return utils.testOverlap(self, otherBlob)

    def merge(self, otherBlob):
        """Merge the given blob into the original blob.

        :param otherBlob:
        :type otherBlob: :class:`pdb_eda.ccp4.DensityBlob`
        """
        self.crsList.update(otherBlob.crsList)
        atoms = self.atoms + [atom for atom in otherBlob.atoms if atom not in self.atoms]
        newBlob = DensityBlob.fromCrsList(self.crsList, self.densityMatrix)

        self.__dict__.update(newBlob.__dict__)
        self.atoms = atoms

    def clone(self):
        """Returns a copy of the density blob.

        :return: densityBlob
        :rtype: :class:`pdb_eda.ccp4.DensityBlob`
        """
        return DensityBlob(self.centroid,self.coordCenter,self.totalDensity,self.volume,self.crsList,self.densityMatrix,self.atoms.copy())

