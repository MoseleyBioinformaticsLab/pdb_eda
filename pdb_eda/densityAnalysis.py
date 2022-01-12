"""
PDB Electron Density Analysis (pdb_eda.densityAnalysis)
-------------------------------------------------------

This module provides methods for the creation of the :class:`pdb_eda.densityAnalysis` class given a PDB id,
along with methods to analyze its 2Fo-Fc and Fo-Fc electron density maps.
"""

import copy
import urllib.request
import os.path
import gzip

import itertools
import collections
import json
import numpy as np
import Bio.PDB as biopdb
import scipy.spatial
from scipy import stats


from . import ccp4
from . import pdbParser

try:
    from . import cutils as utils
except ImportError:
    from . import utils

## Starting data originally from https://arxiv.org/pdf/0804.2488.pdf
paramsPath = os.path.join(os.path.dirname(__file__), 'conf/optimized_params.json')
f000ParamsPath = os.path.join(os.path.dirname(__file__), 'conf/f000_parameters.json.gz')

paramsGlobal = None
with open(paramsPath, 'r') as fh:
    paramsGlobal = json.load(fh)

radiiGlobal = paramsGlobal['radii']
slopesGlobal = paramsGlobal['slopes']
bondedAtomsGlobal = paramsGlobal['bonded_atoms']
elementElectronsGlobal = None
masterFullAtomNameMapElectronsGlobal = None
fullAtomNameMapElectronsGlobal = paramsGlobal['full_atom_name_map_electrons']
fullAtomNameMapAtomTypeGlobal = paramsGlobal['full_atom_name_map_atom_type']
atomTypeLengthGlobal = max(len(atomType) for atomType in fullAtomNameMapAtomTypeGlobal.values()) + 5

def setGlobals(params):
    """Sets global parameters.  Typically used for optimizing parameters.

    :param  params: dictionary of parameters needed for pdb_eda calculations.
    :type params: :py:class:`dict`
    """
    global paramsGlobal
    global radiiGlobal
    global slopesGlobal
    global bondedAtomsGlobal
    global fullAtomNameMapElectronsGlobal
    global fullAtomNameMapAtomTypeGlobal
    global atomTypeLengthGlobal

    paramsGlobal = params
    radiiGlobal = paramsGlobal['radii']
    slopesGlobal = paramsGlobal['slopes']
    bondedAtomsGlobal = paramsGlobal['bonded_atoms']
    fullAtomNameMapElectronsGlobal = paramsGlobal['full_atom_name_map_electrons']
    fullAtomNameMapAtomTypeGlobal = paramsGlobal['full_atom_name_map_atom_type']
    atomTypeLengthGlobal = max(len(atomType) for atomType in fullAtomNameMapAtomTypeGlobal.values()) + 5

def loadF000Parameters():
    """Loads and assigns global parameters needed for F000 estimation."""
    with gzip.open(f000ParamsPath, 'rt') as gzipFile:
        f000Params = json.load(gzipFile)

    global elementElectronsGlobal
    elementElectronsGlobal = f000Params["element_map_electrons"]
    global masterFullAtomNameMapElectronsGlobal
    masterFullAtomNameMapElectronsGlobal = f000Params["full_atom_name_map_electrons"]

ccp4urlPrefix = "http://www.ebi.ac.uk/pdbe/coordinates/files/"
ccp4folder = './ccp4_data/'
pdbfolder = './pdb_data/'
pdburlPrefix = "http://ftp.rcsb.org/pub/pdb/data/structures/all/pdb/"
mmcifurlPrefix = "http://ftp.rcsb.org/pub/pdb/data/structures/all/mmCIF/"


def fromPDBid(pdbid, ccp4density=True, ccp4diff=True, pdbbio=True, pdbi=True, downloadFile=True, mmcif=False):
    """Creates :class:`pdb_eda.densityAnalysis.DensityAnalysis` object given the PDB id if the id is valid
    and the structure has electron density file available.

    :param pdbid: PDB entry ID.
    :type pdbid: :py:class:`str`
    :param ccp4density: Whether to generate ccp4 density object, defaults to :py:obj:`True`.
    :type ccp4density: :py:class:`bool`
    :param ccp4diff: Whether to generate in default of ccp4 difference density object, defaults to :py:obj:`True`.
    :type ccp4diff: :py:class:`bool`
    :param pdbbio: Whether to generate in default of bio.PDB object, defaults to :py:obj:`True`.
    :type pdbbio: :py:class:`bool`
    :param pdbi: Whether to generate in default of PDB object, defaults to :py:obj:`True`.
    :type pdbi: :py:class:`bool`
    :param downloadFile: Whether to save the downloaded ccp4 density, ccp4 difference density, and PDB file, defaults to :py:obj:`True`.
    :type downloadFile: :py:class:`bool`
    :param mmcif: Whether to download the mmCif file, defaults to :py:obj:`False`.
    :type mmcif: :py:class:`bool`

    :return: densityAnalysis object
    :rtype: :class:`pdb_eda.densityAnalysis.DensityAnalysis`
    """
    pdbid = pdbid.lower()

    densityObj = None
    diffDensityObj = None
    pdbObj = None
    biopdbObj = None

    try:
        if ccp4density:
            ## ccp4 2Fo - Fc map parser
            if downloadFile:
                if not os.path.exists(ccp4folder):
                    os.makedirs(ccp4folder)

                ccp4file = ccp4folder + pdbid + '.ccp4'
                if not os.path.isfile(ccp4file):
                    url = ccp4urlPrefix + pdbid + '.ccp4'
                    urllib.request.urlretrieve(url, ccp4file)
                densityObj = ccp4.read(ccp4file, pdbid)
            else:
                densityObj = ccp4.readFromPDBID(pdbid)
            densityObj.densityCutoff = densityObj.meanDensity + 1.5 * densityObj.stdDensity
            densityObj.densityCutoffFromHeader = densityObj.header.densityMean + 1.5 * densityObj.header.rmsd

        if ccp4diff:
            ## ccp4 Fo - Fc map parser
            if downloadFile:
                if not os.path.exists(ccp4folder):
                    os.makedirs(ccp4folder)

                ccp4diffFile = ccp4folder + pdbid + '_diff.ccp4'
                if not os.path.isfile(ccp4diffFile):
                    url = ccp4urlPrefix + pdbid + '_diff.ccp4'
                    urllib.request.urlretrieve(url, ccp4diffFile)

                diffDensityObj = ccp4.read(ccp4diffFile, pdbid)
            else:
                diffDensityObj = ccp4.readFromPDBID(pdbid + '_diff')
            diffDensityObj.diffDensityCutoff = diffDensityObj.meanDensity + 3 * diffDensityObj.stdDensity

        if pdbbio or pdbi:
            pdbfile = pdbfolder + 'pdb' + pdbid + '.ent.gz'
            if not os.path.isfile(pdbfile):
                if not os.path.exists(pdbfolder):
                    os.makedirs(pdbfolder)

                url = pdburlPrefix + 'pdb' + pdbid + '.ent.gz'
                urllib.request.urlretrieve(url, pdbfile)

            if pdbbio:
                # Bio Python PDB parser
                with gzip.open(pdbfile, 'rt') as gzipFile:
                    parser = biopdb.PDBParser(QUIET=True)
                    biopdbObj = parser.get_structure(pdbid, gzipFile)
            if pdbi:
                with gzip.open(pdbfile, 'rt') as gzipFile:
                    pdbObj = pdbParser.readPDBfile(gzipFile)

        if mmcif and downloadFile:
            mmcifFile = pdbfolder + pdbid + '.cif.gz'
            if not os.path.isfile(mmcifFile):
                if not os.path.exists(pdbfolder):
                    os.makedirs(pdbfolder)

                url = mmcifurlPrefix + pdbid + '.cif.gz'
                urllib.request.urlretrieve(url, mmcifFile)
    except:
        return 0

    return DensityAnalysis(pdbid, densityObj, diffDensityObj, biopdbObj, pdbObj)


def fromFile(pdbFile, ccp4DensityFile=None, ccp4DiffDensityFile=None):
    """Creates :class:`pdb_eda.densityAnalysis.DensityAnalysis` object given the appropriate PDB and CCP4 files.

    :param pdbFile: PDB entry file.
    :type pdbFile: :py:class:`str`, :class:`io.IOBase`
    :param ccp4DensityFile: ccp4 density file.
    :type ccp4DensityFile: :py:class:`str`, :class:`io.IOBase`, optional
    :param ccp4DiffDensityFile: ccp4 difference density file.
    :type ccp4DiffDensityFile: :py:class:`str`, :class:`io.IOBase`, optional

    :return: densityAnalysis object
    :rtype: :class:`pdb_eda.densityAnalysis.DensityAnalysis`
    """
    pdbid = "xxxx"
    densityObj = None
    diffDensityObj = None

    try:
        if ccp4DensityFile != None:
            if isinstance(ccp4DensityFile,str):
                densityObj = ccp4.read(ccp4DensityFile, pdbid)
            else:
                densityObj = ccp4.parse(ccp4DensityFile, pdbid)

            densityObj.densityCutoff = densityObj.meanDensity + 1.5 * densityObj.stdDensity
            densityObj.densityCutoffFromHeader = densityObj.header.densityMean + 1.5 * densityObj.header.rmsd

        if ccp4DiffDensityFile != None:
            if isinstance(ccp4DiffDensityFile,str):
                diffDensityObj = ccp4.read(ccp4DiffDensityFile, pdbid)
            else:
                diffDensityObj = ccp4.parse(ccp4DiffDensityFile, pdbid)
            diffDensityObj.diffDensityCutoff = diffDensityObj.meanDensity + 3 * diffDensityObj.stdDensity

        if isinstance(pdbFile,str) and pdbFile.endswith(".gz"):
            with gzip.open(pdbFile, 'rt') as gzipFile:
                parser = biopdb.PDBParser(QUIET=True)
                biopdbObj = parser.get_structure(pdbid, gzipFile)
            with gzip.open(pdbFile, 'rt') as gzipFile:
                pdbObj = pdbParser.readPDBfile(gzipFile)
        else:
            parser = biopdb.PDBParser(QUIET=True)
            biopdbObj = parser.get_structure(pdbid, pdbFile)
            pdbObj = pdbParser.readPDBfile(pdbFile)
    except:
        return 0

    return DensityAnalysis(pdbid, densityObj, diffDensityObj, biopdbObj, pdbObj)


def cleanPDBid(pdbid):
    """Removes PDB entry, ccp4, and mmcif files associated with a PDB id.

    :param pdbid: PDB entry ID.
    :type pdbid: :py:class:`str`

    :return: bool whether the operation was successful.
    :rtype: :py:class:`bool`
    """
    pdbid = pdbid.lower()
    try:
        ccp4file = ccp4folder + pdbid + '.ccp4'
        if os.path.isfile(ccp4file):
            os.remove(ccp4file)

        ccp4diffFile = ccp4folder + pdbid + '_diff.ccp4'
        if os.path.isfile(ccp4diffFile):
            os.remove(ccp4diffFile)

        pdbfile = pdbfolder + 'pdb' + pdbid + '.ent.gz'
        if os.path.isfile(pdbfile):
            os.remove(pdbfile)

        mmcifFile = pdbfolder + pdbid + '.cif.gz'
        if os.path.isfile(mmcifFile):
            os.remove(mmcifFile)
    except:
        return False
    return True

def testCCP4URL(pdbid):
    """Test whether the pdbid has electron density maps by querying if the PDBe API has electron density statistics.
    :param pdbid: PDB entry ID.
    :type pdbid: :py:class:`str`

    :return: bool on test success.
    :rtype: :py:class:`bool`
    """
    try:
        url = "https://www.ebi.ac.uk/pdbe/api/pdb/entry/electron_density_statistics/" + pdbid
        request = urllib.request.urlopen(url)
    except urllib.request.HTTPError:
        return False
    return True


class DensityAnalysis(object):
    """DensityAnalysis class that stores the density, difference density, bio.PDB, and PDB objects."""

    def __init__(self, pdbid, densityObj=None, diffDensityObj=None, biopdbObj=None, pdbObj=None):
        """`densityAnalysis` initializer. Leave `densityObj`, `diffDensityObj`, `biopdbObj` and `pdbObj` as :py:obj:`None`
        to be created. They are not required for initialization but could be required for some methods.

        :param pdbid: PDB entry ID.
        :type pdbid: :py:class:`str`
        :param densityObj: DensityMatrix object.
        :type densityObj: :class:`pdb_eda.ccp4.DensityMatrix`, optional
        :param diffDensityObj: optional DensityMatrix object.
        :type diffDensityObj: :class:`pdb_eda.ccp4.DensityMatrix`, optional
        :param biopdbObj: optional Bio.PDB Structure object.
        :type biopdbObj: :class:`Bio.PDB.Structure.Structure`, optional
        :param pdbObj: optional PDBentry object.
        :type pdbObj: :class:`pdb_eda.pdbParser.PDBentry`, optional
        """
        self.pdbid = pdbid
        self.densityObj = densityObj
        self.diffDensityObj = diffDensityObj
        self.biopdbObj = biopdbObj
        self.pdbObj = pdbObj

        self._symmetryAtoms = None
        self._symmetryOnlyAtoms = None
        self._asymmetryAtoms = None
        self._symmetryAtomCoords = None
        self._symmetryOnlyAtomCoords = None
        self._asymmetryAtomCoords = None
        self._greenBlobList = None
        self._redBlobList = None
        self._blueBlobList = None
        self._fc = None

        self._medians = None
        self._atomCloudDescriptions = None
        self._residueCloudDescriptions = None
        self._domainCloudDescriptions = None
        self._F000 = None
        self._densityElectronRatio = None
        self._numVoxelsAggregated = None
        self._totalAggregatedElectrons = None
        self._totalAggregatedDensity = None

        self._atomTypeOverlapCompleteness = None
        self._atomTypeOverlapIncompleteness = None

    @property
    def symmetryAtoms(self):
        """Returns list of symmetry atoms.

        :return: symmetryAtoms
        :rtype: :py:class:`list`
        """
        if self._symmetryAtoms is None:
            self._calculateSymmetryAtoms()
        return self._symmetryAtoms

    @property
    def symmetryOnlyAtoms(self):
        """Returns list of non-[0,0,0,0] symmetry atoms.

        :return: symmetryAtoms
        :rtype: :py:class:`list`
        """
        if self._symmetryOnlyAtoms is None:
            self._calculateSymmetryAtoms()
        return self._symmetryOnlyAtoms

    @property
    def asymmetryAtoms(self):
        """Returns list of [0,0,0,0] symmetry atoms.

        :return: symmetryAtoms
        :rtype: :py:class:`list`
        """
        if self._asymmetryAtoms is None:
            self._calculateSymmetryAtoms()
        return self._asymmetryAtoms

    @property
    def symmetryAtomCoords(self):
        """Returns list of symmetry atom xyz coordinates.

        :return: symmetryAtomCoords
        :rtype: :py:class:`list`
        """
        if self._symmetryAtoms is None:
            self._calculateSymmetryAtoms()
        return self._symmetryAtomCoords

    @property
    def symmetryOnlyAtomCoords(self):
        """Returns list of non-[0,0,0,0] symmetry atom xyz coordinates.

        :return: symmetryAtomCoords
        :rtype: :py:class:`list`
        """
        if self._symmetryOnlyAtoms is None:
            self._calculateSymmetryAtoms()
        return self._symmetryOnlyAtomCoords

    @property
    def asymmetryAtomCoords(self):
        """Returns list of [0,0,0,0] symmetry atom xyz coordinates.

        :return: symmetryAtomCoords
        :rtype: :py:class:`list`
        """
        if self._asymmetryAtoms is None:
            self._calculateSymmetryAtoms()
        return self._asymmetryAtomCoords

    @property
    def greenBlobList(self):
        """Returns list of all positive significant difference density blobs.

        :return: greenBlobList
        :rtype: :py:class:`list`
        """
        if self._greenBlobList is None:
            self._greenBlobList = self.diffDensityObj.createFullBlobList(self.diffDensityObj.diffDensityCutoff)
        return self._greenBlobList

    @property
    def redBlobList(self):
        """Returns list of all negative significant difference density blobs.

        :return: redBlobList
        :rtype: :py:class:`list`
        """
        if self._redBlobList is None:
            self._redBlobList = self.diffDensityObj.createFullBlobList(-1 * self.diffDensityObj.diffDensityCutoff)
        return self._redBlobList

    @property
    def blueBlobList(self):
        """Returns list of all positive significant density blobs.

        :return: blueBlobList
        :rtype: :py:class:`list`
        """
        if self._blueBlobList is None:
            self._blueBlobList = self.densityObj.createFullBlobList(self.densityObj.densityCutoff)
        return self._blueBlobList

    @property
    def fc(self):
        """Returns the Fc map = 2Fo-Fc - Fo-Fc.

        :return: fc
        :rtype: :class:`pdb_eda.ccp4.DensityMatrix`
        """
        if self._fc is None:
            self._fc = copy.deepcopy(self.densityObj)
            self._fc.density = self.densityObj.density - self.diffDensityObj.density * 2
        return self._fc

    @property
    def fo(self):
        """Returns the Fo map = 2Fo-Fc.

        :return: fo
        :rtype: :class:`pdb_eda.ccp4.DensityMatrix`
        """
        return self.densityObj # using the 2Fo-Fc as the Fo map.

    @property
    def F000(self):
        """Returns estimated F000.  This estimate may not be that accurate.

        :return: F000
        :rtype: :py:class:`float`
        """
        if self._F000 is None:
            self._F000 = self.estimateF000()
        return self._F000

    @property
    def medians(self):
        """Returns median field values calculated per atom type.

        :return: medians
        :rtype: :class:`numpy.array`
        """
        if self._medians is None:
            self.aggregateCloud()
        return self._medians

    @property
    def atomCloudDescriptions(self):
        """Returns aggregated atom cloud descriptions.

        :return: atomCloudDescriptions
        :rtype: :class:`numpy.array`
        """
        if self._atomCloudDescriptions is None:
            self.aggregateCloud()
        return self._atomCloudDescriptions

    @property
    def residueCloudDescriptions(self):
        """Returns aggregated residue cloud descriptions.

        :return: residueCloudDescriptions
        :rtype: :py:class:`list`
        """
        if self._residueCloudDescriptions is None:
            self.aggregateCloud()
        return self._residueCloudDescriptions

    @property
    def domainCloudDescriptions(self):
        """Returns aggregated domain cloud descriptions.

        :return: domainCloudDescriptions
        :rtype: :py:class:`list`
        """
        if self._domainCloudDescriptions is None:
            self.aggregateCloud()
        return self._domainCloudDescriptions

    @property
    def numVoxelsAggregated(self):
        """Returns number of aggregated voxels in cloud analysis.

        :return: numVoxelsAggregated
        :rtype: :py:class:`int`
        """
        if self._numVoxelsAggregated is None:
            self.aggregateCloud()
        return self._numVoxelsAggregated

    @property
    def totalAggregatedElectrons(self):
        """Returns total amount of aggregated electrons in cloud analysis.

        :return: totalAggregatedElectrons
        :rtype: :py:class:`float`
        """
        if self._totalAggregatedElectrons is None:
            self.aggregateCloud()
        return self._totalAggregatedElectrons

    @property
    def totalAggregatedDensity(self):
        """Returns total amount of aggregated density in cloud analysis.

        :return: totalAggregatedDensity
        :rtype: :py:class:`float`
        """
        if self._totalAggregatedDensity is None:
            self.aggregateCloud()
        return self._totalAggregatedDensity

    @property
    def densityElectronRatio(self):
        """Returns the density-electron ratio estimated from cloud analysis.

        :return: densityElectronRatio
        :rtype: :py:class:`float`
        """
        if self._densityElectronRatio == None:
            self.aggregateCloud()
        return self._densityElectronRatio

    @property
    def atomTypeOverlapCompleteness(self):
        """Returns atom-type overlap completeness counts.

        :return: atomTypeOverlapCompleteness
        :rtype: :py:class:`dict`
        """
        if self._atomTypeOverlapCompleteness is None:
            self.aggregateCloud()
        return self._atomTypeOverlapCompleteness

    @property
    def atomTypeOverlapIncompleteness(self):
        """Returns atom-type overlap incompleteness counts.

        :return: atomTypeOverlapIncompleteness
        :rtype: :py:class:`dict`
        """
        if self._atomTypeOverlapIncompleteness is None:
            self.aggregateCloud()
        return self._atomTypeOverlapIncompleteness



    residueCloudHeader = ['chain', 'residue_number', 'residue_name', 'local_density_electron_ratio', 'num_voxels', 'electrons', 'volume', 'centroid_xyz']
    domainCloudHeader = residueCloudHeader
    def aggregateCloud(self, minCloudElectrons=25.0, minTotalElectrons=400.0):
        """Aggregate the electron density map clouds by atom, residue, and domain.
        Calculate and populate `densityAnalysis.densityElectronRatio` and `densityAnalysis.medians` data members.

        :param minCloudElectrons: minimum number of electrons needed to aggregate a residue or domain cloud., defaults to 25.0
        :type minCloudElectrons: :py:class:`float`
        :param minTotalElectrons: mininum number of electrons required to calculate a density-electron ratio., defaults to 400.0
        :type minTotalElectrons: :py:class:`float`
        """
        densityObj = self.densityObj
        biopdbObj = self.biopdbObj

        domainClouds = []
        domainPool = []
        domainList = []
        residueList = []
        atomList = []

        currentRadii = radiiGlobal
        currentSlopes = slopesGlobal
        completelyOverlappedAtomTypes = collections.defaultdict(int)
        incompletelyOverlappedAtomTypes = collections.defaultdict(int)

        allAtomClouds = {}
        centroidDistances = []
        for residue in biopdbObj.get_residues():
            if residue.id[0] != ' ': # skip HETATOM residues.
                continue

            for atom in residue.child_list:
                resAtom = residueAtomName(atom)
                if resAtom not in fullAtomNameMapAtomTypeGlobal.keys() or atom.get_occupancy() == 0:
                    continue

                atomClouds = densityObj.findAberrantBlobs(atom.coord, currentRadii[fullAtomNameMapAtomTypeGlobal[resAtom]], densityObj.densityCutoff)
                allAtomClouds[tuple(atom.coord)] = atomClouds
                if len(atomClouds) > 0:
                    centroidDistances.append(min(np.linalg.norm(atom.coord - i.centroid) for i in atomClouds))
        centroidDistanceCutoff = np.nanmedian(centroidDistances) + 2.5 * np.nanstd(centroidDistances) # ~99% cutoff, but this is calculated across all atom-types.

        for residue in biopdbObj.get_residues():
            if residue.id[0] != ' ': # skip HETATOM residues.
                continue

            residuePool = []
            atomCloudIndeces = {}
            for atom in residue.child_list:
                resAtom = residueAtomName(atom)
                if resAtom not in fullAtomNameMapAtomTypeGlobal.keys() or atom.get_occupancy() == 0:
                    continue

                ## Calculate atom clouds
                atomClouds = allAtomClouds[tuple(atom.coord)]
                if len(atomClouds) == 0:
                    continue
                elif len(atomClouds) == 1:
                    bestAtomCloud = atomClouds[0]
                else:
                    distances = [np.linalg.norm(atom.coord - i.centroid) for i in atomClouds]
                    minDistance = min(distances)
                    if minDistance > centroidDistanceCutoff:
                        continue
                    index = distances.index(minDistance)
                    bestAtomCloud = atomClouds[index]

                for aCloud in atomClouds:
                    aCloud.atoms = [atom]
                atomCloudIndeces[resAtom] = [len(residuePool)+index for index in range(len(atomClouds))]
                residuePool = residuePool + atomClouds ## For aggregating atom clouds into residue clouds

                atomList.append([residue.parent.id, residue.id[1], atom.parent.resname, atom.name, fullAtomNameMapAtomTypeGlobal[resAtom],
                                 bestAtomCloud.totalDensity / fullAtomNameMapElectronsGlobal[resAtom] / atom.get_occupancy(), len(bestAtomCloud.crsList),
                                 fullAtomNameMapElectronsGlobal[resAtom], atom.get_bfactor(), np.linalg.norm(atom.coord - bestAtomCloud.centroid), bestAtomCloud.centroid])
            ## End atom loop

            overlap = np.zeros((len(residuePool), len(residuePool)))
            for i in range(len(residuePool)):
                for j in range(i+1, len(residuePool)):
                    overlap[i][j] = overlap[j][i] = utils.testOverlap(residuePool[i],residuePool[j])


            ## Calculate atom-type overlap completeness.  Needed for parameter optimization.
            for atom in residue.child_list:
                resAtom = residueAtomName(atom)
                if resAtom in atomCloudIndeces:
                    if all(any(overlap[index1][index2] for index1 in atomCloudIndeces[resAtom] for index2 in atomCloudIndeces[resAtom2]) for resAtom2 in bondedAtomsGlobal[resAtom] if resAtom2 in atomCloudIndeces):
                        completelyOverlappedAtomTypes[fullAtomNameMapAtomTypeGlobal[resAtom]] += 1
                    else:
                        incompletelyOverlappedAtomTypes[fullAtomNameMapAtomTypeGlobal[resAtom]] += 1

            ## Group connected residue density clouds together from individual atom clouds
            resClouds = []
            usedIdx = set()
            for startingIndex in range(len(residuePool)):
                if startingIndex not in usedIdx:
                    newCluster = {index for index, o in enumerate(overlap[startingIndex]) if o}
                    currCluster = set([startingIndex])
                    currCluster.update(newCluster)
                    while len(newCluster):
                        newCluster = {index for oldIndex in newCluster for index, o in enumerate(overlap[oldIndex]) if index not in currCluster and o}
                        currCluster.update(newCluster)

                    usedIdx.update(currCluster)
                    resCloud = residuePool[currCluster.pop()].clone()
                    for idx in currCluster:
                        resCloud.merge(residuePool[idx])
                    resClouds.append(resCloud)

            for cloud in resClouds:
                resElectrons = sum([fullAtomNameMapElectronsGlobal[residueAtomName(atom)] * atom.get_occupancy() for atom in cloud.atoms])
                if resElectrons >= minCloudElectrons:
                    residueList.append([residue.parent.id, residue.id[1], residue.resname, cloud.totalDensity / resElectrons, len(cloud.crsList), resElectrons, len(cloud.crsList) * densityObj.header.unitVolume,
                                        cloud.centroid])

            domainPool = domainPool + resClouds ## For aggregating residue clouds into domain clouds
        ## End residue

        ## Group connected domain density clouds together from individual residue clouds
        overlap = np.zeros((len(domainPool), len(domainPool)))
        for i in range(len(domainPool)):
            for j in range(i+1, len(domainPool)):
                overlap[i][j] = overlap[j][i] = utils.testOverlap(domainPool[i],domainPool[j])

        usedIdx = set()
        for startingIndex in range(len(domainPool)):
            if startingIndex not in usedIdx:
                newCluster = {index for index, o in enumerate(overlap[startingIndex]) if o}
                currCluster = set([startingIndex])
                currCluster.update(newCluster)
                while len(newCluster):
                    newCluster = {index for oldIndex in newCluster for index, o in enumerate(overlap[oldIndex]) if index not in currCluster and o}
                    currCluster.update(newCluster)

                usedIdx.update(currCluster)
                domainCloud = domainPool[currCluster.pop()].clone()
                for idx in currCluster:
                    domainCloud.merge(domainPool[idx])
                domainClouds.append(domainCloud)
        ##End domain

        ## Calculate densityElectronRatio, which is technically a weighted mean value now.
        numVoxels = 0
        totalElectrons = 0
        totalDensity = 0
        for cloud in domainClouds:
            atom = cloud.atoms[0]
            domainElectrons = sum([fullAtomNameMapElectronsGlobal[residueAtomName(atom)] * atom.get_occupancy() for atom in cloud.atoms])
            totalElectrons += domainElectrons
            numVoxels += len(cloud.crsList)
            totalDensity += cloud.totalDensity

            if domainElectrons >= minCloudElectrons:
                domainList.append([atom.parent.parent.id, atom.parent.id[1], atom.parent.resname, cloud.totalDensity / domainElectrons, len(cloud.crsList), domainElectrons, len(cloud.crsList) * densityObj.header.unitVolume,
                                  cloud.centroid])

        if totalElectrons < minTotalElectrons:
            return
        else:
            densityElectronRatio = totalDensity / totalElectrons
            domainList.sort(key=lambda x: x[3])
        ## End calculate densityElectronRatio


        def calcSlope(data, atom_type):
            if len(data['chain']) <= 2 or len(np.unique(data['bfactor'])) == 1: ## Less than three data points or all b factors are the same. Should change this to something more robust like 15 unique values.
                return currentSlopes[atom_type]

            slope, intercept, r_vanue, p_value, std_err = stats.linregress(np.log(data['bfactor']), (data['adj_density_electron_ratio']-densityElectronRatio)/densityElectronRatio)
            return currentSlopes[atom_type] if p_value > 0.05 else slope

        try:
            dataType = np.dtype([('chain', np.dtype(('U', 20))), ('residue_number', int), ('residue_name', np.dtype(('U', 10)) ), ('atom_name', np.dtype(('U', 10)) ), ('atom_type', np.dtype(('U', atomTypeLengthGlobal)) ),
                                 ('density_electron_ratio', float), ('num_voxels', int), ('electrons', int), ('bfactor', float), ('centroid_distance', float), ('centroid_xyz', float, (3,)), ('adj_density_electron_ratio', float),
                                 ('domain_fraction', float), ('corrected_fraction', float), ('corrected_density_electron_ratio', float), ('volume', float)])
            atoms = np.asarray([tuple(atom+[0.0 for x in range(5)]) for atom in atomList],dataType) # must pass in list of tuples to create ndarray correctly.
            if not np.isnan(atoms['centroid_distance']).all():
                centroidCutoff = np.nanmedian(atoms['centroid_distance']) + np.nanstd(atoms['centroid_distance']) * 2
                atoms = atoms[atoms['centroid_distance'] < centroidCutoff]
            atom_types = np.unique(atoms['atom_type'])
            medians = { column : {atom_type : np.nanmedian(atoms[column][atoms['atom_type'] == atom_type]) for atom_type in atom_types} for column in ['num_voxels'] }
            medians_translator = np.vectorize(lambda column,atom_type : medians[column][atom_type])

            ## Normalize by volume (numVoxels)
            atoms['adj_density_electron_ratio'] = atoms['density_electron_ratio'] / atoms['num_voxels'] * medians_translator('num_voxels', atoms['atom_type'])
            atoms['volume'] = atoms['num_voxels'] * densityObj.header.unitVolume
            medians.update({column : {atom_type : np.nanmedian(atoms[column][atoms['atom_type'] == atom_type]) for atom_type in atom_types} for column in
                         ['density_electron_ratio', 'centroid_distance', 'adj_density_electron_ratio', 'volume']})
            medians['bfactor'] = {atom_type : np.nanmedian(atoms['bfactor'][(atoms['atom_type'] == atom_type) & (atoms['bfactor'] > 0)]) for atom_type in atom_types}
            atoms['bfactor'][atoms['bfactor'] <= 0] = medians_translator('bfactor', atoms['atom_type'])[atoms['bfactor'] <= 0]
            medians['slopes'] = {atom_type : calcSlope(atoms[atoms['atom_type'] == atom_type], atom_type) for atom_type in atom_types}

            ## Correct by b-factor
            atoms['domain_fraction'] = (atoms['adj_density_electron_ratio'] - densityElectronRatio) / densityElectronRatio
            atoms['corrected_fraction'] = atoms['domain_fraction'] - (np.log(atoms['bfactor']) - np.log(medians_translator('bfactor', atoms['atom_type']))) * medians_translator('slopes', atoms['atom_type'])
            atoms['corrected_density_electron_ratio'] = atoms['corrected_fraction'] * densityElectronRatio + densityElectronRatio
            medians.update({column : {atom_type : np.nanmedian(atoms[column][atoms['atom_type'] == atom_type]) for atom_type in atom_types} for column in
                         ['domain_fraction', 'corrected_fraction', 'corrected_density_electron_ratio']})
        except:
            return

        self._densityElectronRatio = densityElectronRatio
        self._numVoxelsAggregated = numVoxels
        self._totalAggregatedElectrons = totalElectrons
        self._totalAggregatedDensity = totalDensity
        self._medians = medians
        self._atomCloudDescriptions = atoms
        self._residueCloudDescriptions = residueList
        self._domainCloudDescriptions = domainList
        self._atomTypeOverlapCompleteness = completelyOverlappedAtomTypes
        self._atomTypeOverlapIncompleteness = incompletelyOverlappedAtomTypes


    def medianAbsFoFc(self):
        """Calculates median absolute values for the Fo and Fc maps less than 1 sigma.
        These values should be comparable, i.e. low relative difference, for RSCC and RSR metric calculations.

        :return: tuple of median abs values from fo and fc maps respectively.
        :rtype: :py:class:`tuple`
        """
        fo = self.fo
        fc = self.fc
        foDensityCutoff = fo.meanDensity + 1.0 * fo.stdDensity
        fcDensityCutoff = fc.meanDensity + 1.0 * fc.stdDensity
        ncrs = fo.header.uniqueNcrs
        densityPairs = [ (density, diffDensity) for density, diffDensity in ((utils.getPointDensityFromCrs(fo, crs),utils.getPointDensityFromCrs(fc, crs)) for crs in itertools.product(range(ncrs[0]),range(ncrs[1]),range(ncrs[2])))
          if abs(density) < foDensityCutoff and abs(diffDensity) < fcDensityCutoff ]

        foValues = np.asarray([pair[0] for pair in densityPairs])
        fcValues = np.asarray([pair[1] for pair in densityPairs])
        return (np.median(np.abs(foValues)),np.median(np.abs(fcValues)))

    residueMetricsHeaderList = ['chain', 'residue_number', 'residue_name', "rscc", "rsr", "mean_occupancy", "occupancy_weighted_mean_bfactor"]
    def residueMetrics(self, residueList=None):
        """RETURNS rscc and rsr statistics for each residue using the Fo and Fc density maps.

        :param residueList:
        :type residueList: :py:class:`list`, optional

        :return: results
        :rtype: :py:class:`list`
        """
        resolution = self.biopdbObj.header['resolution']
        radius = 0.7
        if 0.6 <= resolution <= 3:
            radius = (resolution - 0.6) / 3 + 0.7
        elif resolution > 3:
            radius = resolution * 0.5

        if residueList == None:
            residueList = list(self.biopdbObj.get_residues())

        residueResults = []
        for residue in residueList:
            crsList = set()
            bfactorWeightedSum = occupancySum = 0.0
            for atom in residue.child_list:
                crsList.update(utils.getSphereCrsFromXyz(self.fo, atom.coord, radius, 0.0))
                bfactorWeightedSum += atom.get_bfactor() * atom.get_occupancy()
                occupancySum += atom.get_occupancy()

            (rscc, rsr) = self.calculateRsccRsrMetrics(crsList)
            residueResults.append([residue.parent.id, residue.id[1], residue.resname, rscc, rsr, occupancySum / len(residue.child_list), bfactorWeightedSum / occupancySum])

        return residueResults

    atomMetricsHeaderList = ['chain', 'residue_number', 'residue_name', "atom_name", "symmetry", "xyz", "rscc", "rsr", "occupancy", "bfactor"]
    def atomMetrics(self, atomList=None):
        """RETURNS rscc and rsr statistics for each residue using the Fo and Fc density maps.

        :param atomList:
        :type atomList: :py:class:`list`, optional

        :return: results
        :rtype: :py:class:`list`
        """
        resolution = self.biopdbObj.header['resolution']
        radius = 0.7
        if 0.6 <= resolution <= 3:
            radius = (resolution - 0.6) / 3 + 0.7
        elif resolution > 3:
            radius = resolution * 0.5

        if atomList == None:
            atomList = self.asymmetryAtoms

        atomResults = []
        for atom in atomList:
            crsList = set(utils.getSphereCrsFromXyz(self.fo, atom.coord, radius, 0.0))
            (rscc, rsr) = self.calculateRsccRsrMetrics(crsList)
            atomResults.append([atom.parent.parent.id, atom.parent.id[1], atom.parent.resname, atom.name, atom.symmetry, atom.coord, rscc, rsr, atom.get_occupancy(), atom.get_bfactor()])

        return atomResults

    def calculateRsccRsrMetrics(self, crsList):
        """Calculates and returns RSCC and RSR metrics.
        This method of calculating RSCC and RSR assumes that the Fo and Fc maps are appropriately scaled.
        Comparison of median absolute values below one sigma should be quite similar between Fo and Fc maps.

        :param crsList:
        :type crsList: :py:class:`list`, :py:class:`set`

        :return: rscc_rsr_arrays_tuple
        :rtype: :py:obj:`tuple`
        """
        fo = self.fo
        fc = self.fc
        foDensity = np.asarray([utils.getPointDensityFromCrs(fo,i) for i in crsList])
        fcDensity = np.asarray([utils.getPointDensityFromCrs(fc, i) for i in crsList])

        rscc = stats.stats.pearsonr(foDensity, fcDensity)[0]
        rsr = sum(abs(foDensity - fcDensity)) / sum(abs(foDensity + fcDensity))
        return (rscc,rsr)


    def _calculateSymmetryAtoms(self):
        """Calculate all the symmetry and nearby cells and keep those have at least on atom within 5 grid points of the non-repeating crs boundary.
        Ref: Biomolecular Crystallography: Principles, Practice, and Application to Structural Biology by Bernhard Rupp.
        Orthogonalization matrix O and deororthogonalization matrix O' are from :class:`pdb_eda.ccp4` object.
        Rotation matrix R and Translation matrix T is from :class:`pdb_eda.pdbParser` object.
        The neighbouring cells can be calculated using formula,
        X' = O(O'(RX + T) + T') = OO'(RX+T) + OT' = RX + T + O[-1/0/1,-1/0/1,-1/0/1].
        Assign the list of :class:`pdb_eda.densityAnalysis.symAtom` instances to `densityAnalysis.symmetryAtoms` data member
        """
        densityObj = self.densityObj
        biopdbObj = self.biopdbObj
        pdbObj = self.pdbObj

        ## For inRangeAtoms, the min/max range of xyz axes (the circumscribed box)
        ncrs = densityObj.header.ncrs
        orginalDensityBox = [densityObj.header.crs2xyzCoord(i) for i in [[c, r, s] for c in [0, ncrs[0]-1] for r in [0, ncrs[1]-1] for s in [0, ncrs[2]-1]]]
        xs = sorted([i[0] for i in orginalDensityBox])
        ys = sorted([i[1] for i in orginalDensityBox])
        zs = sorted([i[2] for i in orginalDensityBox])

        allAtoms = utils.createSymmetryAtoms(list(biopdbObj.get_atoms()), pdbObj.header.rotationMats, densityObj.header.orthoMat, xs,ys,zs)

        self._symmetryAtoms = allAtoms
        self._symmetryAtomCoords = np.asarray([atom.coord for atom in allAtoms])
        self._symmetryOnlyAtoms = [atom for atom in allAtoms if atom.symmetry != (0,0,0,0)]
        self._symmetryOnlyAtomCoords = np.asarray([atom.coord for atom in self._symmetryOnlyAtoms])
        self._asymmetryAtoms = [atom for atom in allAtoms if atom.symmetry == (0,0,0,0)]
        self._asymmetryAtomCoords = np.asarray([atom.coord for atom in self._asymmetryAtoms])

    blobStatisticsHeader = ['distance_to_atom', 'sign', 'electrons_of_discrepancy', 'num_voxels', 'volume', 'chain', 'residue_number', 'residue_name', 'atom_name', 'atom_symmetry', 'atom_xyz', 'centroid_xyz']
    def calculateAtomSpecificBlobStatistics(self, blobList):
        """Calculate atom-specific blob statistics.

        :param blobList: list of blobs to calculate statistics for.
        :type blobList: :py:class:`list`

        :return blobStats: Difference density map statistics.
        :rtype: :py:class:`list`
        """
        symmetryAtoms = self.symmetryAtoms
        symmetryAtomCoords = self.symmetryAtomCoords

        if not self.densityElectronRatio:
            raise RuntimeError("Failed to calculate densityElectronRatio, probably due to total aggregated electrons less than the minimum.")
        densityElectronRatio = self.densityElectronRatio

        blobStats = []
        for blob in blobList:
            centroid = np.array(blob.centroid).reshape((1, 3))
            symmetryDistances = scipy.spatial.distance.cdist(centroid, symmetryAtomCoords)
            atom = symmetryAtoms[np.argmin(symmetryDistances[0])] # closest atom
            sign = '+' if blob.totalDensity >= 0 else '-'
            blobStats.append([symmetryDistances.min(), sign, abs(blob.totalDensity / densityElectronRatio), len(blob.crsList), blob.volume, atom.parent.parent.id, atom.parent.id[1], atom.parent.resname, atom.name, atom.symmetry, atom.coord, blob.centroid])

        return blobStats

    # Headers that match the order of the results
    regionDensityHeader = [ "actual_significant_regional_density", "num_electrons_actual_significant_regional_density" ]
    atomRegionDensityHeader = ['model', 'chain', 'residue_number', 'residue_name', "atom_name", "occupancy"] + regionDensityHeader
    symmetryAtomRegionDensityHeader = ['model', 'chain', 'residue_number', 'residue_name', "atom_name", "symmetry", "atom_xyz", "fully_within_density_map"] + regionDensityHeader
    residueRegionDensityHeader = ['model', 'chain', 'residue_number', 'residue_name', "mean_occupancy"] + regionDensityHeader


    def calculateAtomRegionDensity(self, radius, numSD=1.5, type="", useOptimizedRadii=False):
        """Calculates significant region density in a given radius of each atom.

        :param radius: the search radius.
        :type radius: :py:class:`float`
        :param numSD: number of standard deviations of significance, defaults to 1.5
        :type numSD: py:class:`float`
        :param type: atom type to filter on., defaults to ""
        :type type: :py:class:`str`

        :return diffMapRegionStats: density map region statistics per atom.
        :rtype: :py:class:`list`
        """
        biopdbObj = self.biopdbObj
        atoms = list(biopdbObj.get_atoms())
        if type:
            atoms = [atom for atom in atoms if atom.name == type]

        results = []
        for atom in atoms:
            resAtom = residueAtomName(atom)
            testRadius = radiiGlobal[fullAtomNameMapAtomTypeGlobal[resAtom]] if useOptimizedRadii and resAtom in fullAtomNameMapAtomTypeGlobal.keys() else radius
            result = self.calculateRegionDensity([atom.coord], testRadius, numSD)
            results.append([atom.parent.parent.parent.id, atom.parent.parent.id, atom.parent.id[1], atom.parent.resname, atom.name, atom.get_occupancy()] + result)

        return results

    def calculateSymmetryAtomRegionDensity(self, radius, numSD=1.5, type="", useOptimizedRadii=False):
        """Calculates significant region density in a given radius of each symmetry atom.

        :param radius: the search radius.
        :type radius: :py:class:`float`
        :param numSD: number of standard deviations of significance., defaults to 1.5
        :type numSD: :py:class:`float`
        :param type: atom type to filter on., defaults to ""
        :type type: :py:class:`str`

        :return diffMapRegionStats: density map region statistics per atom.
        :rtype: :py:class:`list`
        """
        atoms = self.symmetryAtoms
        if type:
            atoms = [atom for atom in atoms if atom.name == type]

        results = []
        for atom in atoms:
            resAtom = residueAtomName(atom)
            testRadius = radiiGlobal[fullAtomNameMapAtomTypeGlobal[resAtom]] if useOptimizedRadii and resAtom in fullAtomNameMapAtomTypeGlobal.keys() else radius
            (result,valid) = self.calculateRegionDensity([atom.coord], testRadius, numSD, testValidCrs=True)
            results.append([atom.parent.parent.parent.id, atom.parent.parent.id, atom.parent.id[1], atom.parent.resname, atom.name, atom.symmetry, atom.coord, valid] + result)

        return results

    def calculateResidueRegionDensity(self, radius, numSD=1.5, type="", atomMask=None, useOptimizedRadii=False):
        """Calculates significant region density in a given radius of each residue.

        :param radius: the search radius.
        :type radius: :py:class:`float`
        :param numSD: number of standard deviations of significance., defaults to 1.5
        :type numSD: :py:class:`float`
        :param type: atom type to filter on., defaults to ""
        :type type: :py:class:`str`
        :param atomMask: residue specific atom mask.
        :type type: :py:class:`dict`, optional

        :return diffMapRegionStats: density map region statistics per residue.
        :rtype: :py:class:`list`
        """
        biopdbObj = self.biopdbObj

        results = []
        residues = list(biopdbObj.get_residues())
        if type:
            residues = [residue for residue in residues if residue.resname == type]
        for residue in residues:
            atoms = [atom for atom in residue.get_atoms() if not atomMask or residue.resname not in atomMask or atom.name in atomMask[residue.resname]]
            if atoms:
                xyzCoordList = [atom.coord for atom in atoms]
                meanOccupancy = np.mean([atom.get_occupancy() for atom in atoms])
                if useOptimizedRadii:
                    resAtoms = [residueAtomName(atom) for atom in atoms]
                    radii = [radiiGlobal[fullAtomNameMapAtomTypeGlobal[resAtom]] if resAtom in fullAtomNameMapAtomTypeGlobal.keys() else radius for resAtom in resAtoms]
                    result = self.calculateRegionDensity(xyzCoordList, radii, numSD)
                else:
                    result = self.calculateRegionDensity(xyzCoordList, radius, numSD)
                results.append([residue.parent.parent.id, residue.parent.id, residue.id[1], residue.resname, meanOccupancy ] + result)

        return results

    def calculateRegionDensity(self, xyzCoordList, radius, numSD=1.5, testValidCrs=False):
        """Calculate region-specific density from the electron density matrix.

        :param xyzCoordLists: single xyz coordinate or a list of xyz coordinates.
        :type xyzCoordList: :py:class:`list`
        :param radius: the search radius or list of search radii.
        :type radius: :py:class:`float` or :py:class:`list`
        :param numSD: number of standard deviations of significance., defaults to 1.5
        :type numSD: :py:class:`float`
        :param testValidCrs: whether to test crs are valid and return the results., defaults to :py:obj:`False`
        :type testValidCrs: :py:class:`bool`

        :return diffMapRegionStats: density map region statistics and optional validCrs result.
        :rtype: :py:class:`list`, :py:class:`tuple`
        """
        if not self.densityElectronRatio:
            raise RuntimeError("Failed to calculate densityElectronRatio, probably due to total aggregated electrons less than the minimum.")
        densityElectronRatio = self.densityElectronRatio

        densityObj = self.densityObj
        densityCutoff = densityObj.meanDensity + numSD * densityObj.stdDensity

        # observed significant regional density
        blue = densityObj.findAberrantBlobs(xyzCoordList, radius, densityCutoff)
        actual_sig_regional_density = sum([blob.totalDensity for blob in blue])
        num_electrons_actual_sig_regional_density = actual_sig_regional_density / densityElectronRatio

        result = [ actual_sig_regional_density, num_electrons_actual_sig_regional_density ]
        if testValidCrs:
            return (result, utils.testValidXyzList(densityObj, xyzCoordList, radius))
        else:
            return result


    # Headers that match the order of the results
    regionDiscrepancyHeader = [ "actual_abs_significant_regional_discrepancy", "num_electrons_actual_abs_significant_regional_discrepancy",
                 "expected_abs_significant_regional_discrepancy", "num_electrons_expected_abs_significant_regional_discrepancy",
                 "actual_significant_regional_discrepancy", "num_electrons_actual_significant_regional_discrepancy",
                 "actual_positive_significant_regional_discrepancy", "num_electrons_actual_positive_significant_regional_discrepancy",
                 "actual_negative_significant_regional_discrepancy", "num_electrons_actual_negative_significant_regional_discrepancy" ]
    atomRegionDiscrepancyHeader = ['model', 'chain', 'residue_number', 'residue_name', "atom_name", "occupancy"] + regionDiscrepancyHeader
    symmetryAtomRegionDiscrepancyHeader = ['model', 'chain', 'residue_number', 'residue_name', "atom_name", "symmetry", "atom_xyz", "fully_within_density_map"] + regionDiscrepancyHeader
    residueRegionDiscrepancyHeader = ['model', 'chain', 'residue_number', 'residue_name', "mean_occupancy"] + regionDiscrepancyHeader

    def calculateAtomRegionDiscrepancies(self, radius, numSD=3.0, type=""):
        """Calculates significant region discrepancies in a given radius of each atom.

        :param radius: the search radius.
        :type radius: :py:class:`float`
        :param numSD: number of standard deviations of significance, defaults to 3.0
        :type numSD: py:class:`float`
        :param type: atom type to filter on., defaults to ""
        :type type: :py:class:`str`

        :return diffMapRegionStats: Difference density map region statistics per atom.
        :rtype: :py:class:`list`
        """
        biopdbObj = self.biopdbObj
        atoms = list(biopdbObj.get_atoms())
        if type:
            atoms = [atom for atom in atoms if atom.name == type]

        results = []
        for atom in atoms:
            result = self.calculateRegionDiscrepancy([atom.coord], radius, numSD)
            results.append([atom.parent.parent.parent.id, atom.parent.parent.id, atom.parent.id[1], atom.parent.resname, atom.name, atom.get_occupancy()] + result)

        return results

    def calculateSymmetryAtomRegionDiscrepancies(self, radius, numSD=3.0, type=""):
        """Calculates significant region discrepancies in a given radius of each symmetry atom.

        :param radius: the search radius.
        :type radius: :py:class:`float`
        :param numSD: number of standard deviations of significance., defaults to 3.0
        :type numSD: :py:class:`float`
        :param type: atom type to filter on., defaults to ""
        :type type: :py:class:`str`

        :return diffMapRegionStats: Difference density map region statistics per atom.
        :rtype: :py:class:`list`
        """
        atoms = self.symmetryAtoms
        if type:
            atoms = [atom for atom in atoms if atom.name == type]

        results = []
        for atom in atoms:
            (result,valid) = self.calculateRegionDiscrepancy([atom.coord], radius, numSD, testValidCrs=True)
            results.append([atom.parent.parent.parent.id, atom.parent.parent.id, atom.parent.id[1], atom.parent.resname, atom.name, atom.symmetry, atom.coord, valid] + result)

        return results

    def calculateResidueRegionDiscrepancies(self, radius, numSD=3.0, type="", atomMask=None):
        """Calculates significant region discrepancies in a given radius of each residue.

        :param radius: the search radius.
        :type radius: :py:class:`float`
        :param numSD: number of standard deviations of significance., defaults to 3.0
        :type numSD: :py:class:`float`
        :param type: atom type to filter on., defaults to ""
        :type type: :py:class:`str`
        :param atomMask: residue specific atom mask.
        :type type: :py:class:`dict`, optional

        :return diffMapRegionStats: Difference density map region statistics per residue.
        :rtype: :py:class:`list`
        """
        biopdbObj = self.biopdbObj

        results = []
        residues = list(biopdbObj.get_residues())
        if type:
            residues = [residue for residue in residues if residue.resname == type]
        for residue in residues:
            atoms = [atom for atom in residue.get_atoms() if not atomMask or (residue.resname in atomMask and atom.name in atomMask[residue.resname])]
            xyzCoordList = [atom.coord for atom in atoms]
            meanOccupancy = np.mean([atom.get_occupancy() for atom in atoms])
            result = self.calculateRegionDiscrepancy(xyzCoordList, radius, numSD)
            results.append([residue.parent.parent.id, residue.parent.id, residue.id[1], residue.resname, meanOccupancy ] + result)

        return results

    def calculateRegionDiscrepancy(self, xyzCoordList, radius, numSD=3.0, testValidCrs=False):
        """Calculate region-specific discrepancy from the difference density matrix.

        :param xyzCoordLists: single xyz coordinate or a list of xyz coordinates.
        :type xyzCoordList: :py:class:`list`
        :param radius: the search radius.
        :type radius: :py:class:`float`
        :param numSD: number of standard deviations of significance., defaults to 3.0
        :type numSD: :py:class:`float`
        :param testValidCrs: whether to test crs are valid and return the results., defaults to :py:obj:`False`
        :type testValidCrs: :py:class:`bool`

        :return diffMapRegionStats: Difference density map region statistics and optional validCrs result.
        :rtype: :py:class:`list`, :py:class:`tuple`
        """
        if not self.densityElectronRatio:
            raise RuntimeError("Failed to calculate densityElectronRatio, probably due to total aggregated electrons less than the minimum.")
        densityElectronRatio = self.densityElectronRatio

        diffDensityObj = self.diffDensityObj
        diffDensityCutoff = diffDensityObj.meanDensity + numSD * diffDensityObj.stdDensity

        # observed significant regional discrepancy
        green = diffDensityObj.findAberrantBlobs(xyzCoordList, radius, diffDensityCutoff)
        red = diffDensityObj.findAberrantBlobs(xyzCoordList, radius, -1.0 * diffDensityCutoff)
        actual_positive_sig_regional_discrep = sum([blob.totalDensity for blob in green])
        num_electrons_actual_positive_sig_regional_discrep = actual_positive_sig_regional_discrep / densityElectronRatio
        actual_negative_sig_regional_discrep = sum([blob.totalDensity for blob in red])
        num_electrons_actual_negative_sig_regional_discrep = actual_negative_sig_regional_discrep / densityElectronRatio
        actual_sig_regional_discrep = actual_positive_sig_regional_discrep + actual_negative_sig_regional_discrep
        num_electrons_actual_sig_regional_discrep = actual_sig_regional_discrep / densityElectronRatio
        actual_abs_sig_regional_discrep = abs(actual_positive_sig_regional_discrep) + abs(actual_negative_sig_regional_discrep)
        num_electrons_actual_abs_sig_regional_discrep = actual_abs_sig_regional_discrep / densityElectronRatio

        # expected absolute significant regional discrepancy
        total_abs_sig_discrep = diffDensityObj.getTotalAbsDensity(diffDensityCutoff)
        total_voxel_count = len(diffDensityObj.densityArray)
        avg_abs_vox_discrep = total_abs_sig_discrep / total_voxel_count
        regional_voxel_count = len(utils.getSphereCrsFromXyzList(diffDensityObj, xyzCoordList, radius))
        expected_abs_sig_regional_discrep = avg_abs_vox_discrep * regional_voxel_count
        num_electrons_expected_abs_sig_regional_discrep = expected_abs_sig_regional_discrep / densityElectronRatio

        result = [ actual_abs_sig_regional_discrep, num_electrons_actual_abs_sig_regional_discrep,
                 expected_abs_sig_regional_discrep, num_electrons_expected_abs_sig_regional_discrep,
                 actual_sig_regional_discrep, num_electrons_actual_sig_regional_discrep,
                 actual_positive_sig_regional_discrep, num_electrons_actual_positive_sig_regional_discrep,
                 actual_negative_sig_regional_discrep, num_electrons_actual_negative_sig_regional_discrep ]

        if testValidCrs:
            return (result, utils.testValidXyzList(diffDensityObj, xyzCoordList, radius))
        else:
            return result


    def estimateF000(self):
        """Estimate the F000 term as sum of all electrons over the unit cell volume

        :return: estimatedF000
        :rtype: :py:class:`float`
        """
        densityObj = self.densityObj
        biopdbObj = self.biopdbObj
        pdbObj = self.pdbObj

        if not elementElectronsGlobal:
            loadF000Parameters()

        totalElectrons = 0

        allAtoms = list(biopdbObj.get_atoms())
        for atom in allAtoms:
            fullAtomName = residueAtomName(atom)
            if fullAtomName in masterFullAtomNameMapElectronsGlobal:
                totalElectrons += masterFullAtomNameMapElectronsGlobal[fullAtomName]
            elif atom.element in elementElectronsGlobal:
                totalElectrons += elementElectronsGlobal[atom.element] + 1 # + 1 as an estimate for the number of H.

        totalElectrons *= len(pdbObj.header.rotationMats)
        asuVolume = densityObj.header.unitVolume * densityObj.header.nintervalX * densityObj.header.nintervalY * densityObj.header.nintervalZ

        return totalElectrons/asuVolume


def residueAtomName(atom):
    """Returns a combined residue and atom name used to select an atom type.

    :param atom:
    :type atom: :class:`Bio.PDB.atom`

    :return: fullAtomName
    :rtype: :py:class:`str`
    """
    return atom.parent.resname.strip() + '_' + atom.name