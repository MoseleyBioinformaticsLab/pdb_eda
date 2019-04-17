# !/usr/bin/python3

"""
PDB Electron Density Analysis (pdb_eda.densityAnalysis)
-------------------------------------------------------

This module provides methods for the creation of the :class:`pdb_eda.densityAnalysis` class given a PDB id,
along with methods to analyze its electron density.
"""

import copy
import urllib.request
import os.path

import json
import pandas
import numpy as np
import Bio.PDB as biopdb
import scipy.spatial
from scipy import stats

from . import ccp4
from . import pdbParser
from . import validationStats

## Data originally from https://arxiv.org/pdf/0804.2488.pdf
radiiParamPath = os.path.join(os.path.dirname(__file__), 'conf/original_radii_slope_param.json')
electronParamPath = os.path.join(os.path.dirname(__file__), 'conf/atom_type_electron_param.json')

with open(radiiParamPath, 'r') as fh:
    radiiParams = json.load(fh)
with open(electronParamPath, 'r') as fh:
    electronParams = json.load(fh)

radiiDefault = radiiParams['radii']
slopesDefault = radiiParams['slopes']
elementElectrons = electronParams['elementElectrons']
aaElectrons = electronParams['aaElectrons']
electrons = electronParams['atomTypeElectrons']
atomType = electronParams['atomType']

ccp4urlPrefix = "http://www.ebi.ac.uk/pdbe/coordinates/files/"
ccp4urlSuffix = ".ccp4"
ccp4folder = './ccp4_data/'
pdbfolder = './pdb_data/'


def fromPDBid(pdbid, ccp4density=True, ccp4diff=True, pdbbio=True, pdbi=True, downloadFile=True):
    """
    Creates :class:`pdb_eda.densityAnalysis.DensityAnalysis` object given the PDB id if the id is valid
    and the structure has electron density file available.

    :param str pdbid: PDB id.
    :param ccp4density: Whether to generate ccp4 density object. Default is true.
    :param ccp4diff: Whether to generate in default of ccp4 difference density object. Default is true.
    :param pdbbio: Whether to generate in default of bio.PDB object. Default is true.
    :param pdbi: Whether to generate in default of PDB object. Default is true.
    :param downloadFile: Whether to save the downloaded ccp4 density, ccp4 difference density, and PDB file. Default is true.

    :return: :class:`pdb_eda.densityAnalysis`

    :type ccp4density: :py:obj:`True` or :py:obj:`False`
    :type ccp4diff: :py:obj:`True` or :py:obj:`False`
    :type pdbbio: :py:obj:`True` or :py:obj:`False`
    :type pdbi: :py:obj:`True` or :py:obj:`False`
    :type downloadFile: :py:obj:`True` or :py:obj:`False`
    """
    pdbid = pdbid.lower()
    #print("working on " + pdbid + ', ', str(datetime.datetime.now()))

    try:
        if ccp4density:
            ## ccp4 2Fo - Fc map parser
            if downloadFile:
                if not os.path.exists(ccp4folder):
                    os.makedirs(ccp4folder)

                ccp4file = ccp4folder + pdbid + '.ccp4'
                if not os.path.isfile(ccp4file):
                    url = ccp4urlPrefix + pdbid + ccp4urlSuffix
                    urllib.request.urlretrieve(url, ccp4file)
                densityObj = ccp4.read(ccp4file, pdbid)
            else:
                densityObj = ccp4.readFromPDBID(pdbid)
            densityObj.densityCutoff = np.mean(densityObj.densityArray) + 1.5 * np.std(densityObj.densityArray)
            densityObj.densityCutoffFromHeader = densityObj.header.densityMean + 1.5 * densityObj.header.rmsd

            '''
            sample = np.random.choice(densityObj.densityArray, int(len(densityObj.densityArray) / 10))
            kernel = stats.gaussian_kde(sample)
            #kernel = stats.gaussian_kde(densityObj.densityArray)
            x = np.linspace(min(densityObj.densityArray), max(densityObj.densityArray), 200)
            mode = x[np.argmax(kernel(x))]
            leftside = [i for i in densityObj.densityArray if i < mode]
            dev = np.sqrt(sum([(i - mode) ** 2 for i in leftside]) / len(leftside))
            densityObj.densityCutoffFromLeftSide = mode + dev * 1.5
            densityObj.densityCutoffFromLeftSide2 = mode + dev * 2
            densityObj.densityCutoffFromLeftSide25 = mode + dev * 2.5
            densityObj.densityCutoffFromLeftSide3 = mode + dev * 3
            '''

        if ccp4diff:
            ## ccp4 Fo - Fc map parser
            if downloadFile:
                if not os.path.exists(ccp4folder):
                    os.makedirs(ccp4folder)

                ccp4diffFile = ccp4folder + pdbid + '_diff.ccp4'
                if not os.path.isfile(ccp4diffFile):
                    url = ccp4urlPrefix + pdbid + '_diff' + ccp4urlSuffix
                    urllib.request.urlretrieve(url, ccp4diffFile)

                diffDensityObj = ccp4.read(ccp4diffFile, pdbid)
            else:
                diffDensityObj = ccp4.readFromPDBID(pdbid + '_diff')
            diffDensityObj.diffDensityCutoff = np.mean(diffDensityObj.densityArray) + 3 * np.std(diffDensityObj.densityArray)

        if pdbbio or pdbi:
            pdbfile = pdbfolder + 'pdb' + pdbid + '.ent'
            if not os.path.isfile(pdbfile):
                if not os.path.exists(pdbfolder):
                    os.makedirs(pdbfolder)

                pdbl = biopdb.PDBList()
                pdbl.retrieve_pdb_file(pdbid, pdir=pdbfolder, file_format="pdb")

            if pdbbio:
                # Bio Python PDB parser
                parser = biopdb.PDBParser(QUIET=True)
                biopdbObj = parser.get_structure(pdbid, pdbfile)
            if pdbi:
                ## my own PDB parser
                pdbObj = pdbParser.readPDBfile(pdbfile)

        if not downloadFile and os.path.isfile(pdbfile):
            os.remove(pdbfile)
    except:
        return 0

    return DensityAnalysis(pdbid, densityObj, diffDensityObj, biopdbObj, pdbObj)


class DensityAnalysis(object):
    """DensityAnalysis class that stores the density, difference density, bio.PDB, and PDB objects."""

    def __init__(self, pdbid, densityObj=None, diffDensityObj=None, biopdbObj=None, pdbObj=None):
        """
        `densityAnalysis` initializer. Leave `densityObj`, `diffDensityObj`, `biopdbObj` and `pdbObj` as :py:obj:`None`
        to be created. They are not required for initialization but could be required for some methods.

        :param str pdbid: PDB id.
        :param densityObj: Optional :class:`pdb_eda.ccp4` object.
        :param diffDensityObj: Optional :class:`pdb_eda.ccp4` object.
        :param biopdbObj: Optional `bio.PDB` object.
        :param pdbObj: Optional :class:`pdb_eda.pdbParser.PDBentry` object.
        """
        self.pdbid = pdbid
        self.densityObj = densityObj
        self.diffDensityObj = diffDensityObj
        self.biopdbObj = biopdbObj
        self.pdbObj = pdbObj

        self.symmetryAtoms = None
        self.greenBlobList = None
        self.redBlobList = None
        self.chainMedian = None
        self.medians = None
        self.atomList = None
        self.residueList = None
        self.chainList = None
        self.statistics = None
        self.f000 = None
        self.chainNvoxel = None
        self.chainTotalE = None


    def validation(self, densityObj=None, diffDensityObj=None, biopdbObj=None, recalculate=False):
        """
        Populate `DensityAnalysis.statistics` data member with RSR and RSCC.
        Leave `densityObj`, `diffDensityObj`, `biopdbObj` and `pdbObj` as :py:obj:`None` to be read in,
        and it will use its own data member.

        :param str pdbid: PDB id
        :param densityObj: Optional :class:`pdb_eda.ccp4` object.
        :param diffDensityObj: Optional :class:`pdb_eda.ccp4` object.
        :param biopdbObj: Optional `bio.PDB` object.
        :param pdbObj: Optional :class:`pdb_eda.pdbParser.PDBentry` object.
        :param recalculate: Whether or not to recalculate if `densityAnalysis.statistics` already exist.
        :type recalculate: :py:obj:`True` or :py:obj:`False`

        :return: :py:obj:`None`
        """
        if self.statistics and not recalculate:
            return None
        if not densityObj:
            densityObj = self.densityObj
        if not diffDensityObj:
            diffDensityObj = self.diffDensityObj
        if not biopdbObj:
            biopdbObj = self.biopdbObj

        valid = validationStats.validationStats(self.pdbid)
        fo = copy.deepcopy(densityObj)
        fc = copy.deepcopy(densityObj)

        fc.density = densityObj.density - diffDensityObj.density * 2
        sigma3 = 0

        self.statistics = valid.getStats(biopdbObj, fc, fo, sigma3)


    def aggregateCloud(self, radiiUpdate={}, slopesUpdate={}, densityObj=None, biopdbObj=None, atomL=False, residueL=False, chainL=False, recalculate=False):
        """
        Aggregate the electron density map clouds by atom, residue, and chain.
        Calculate and populate `densityAnalysis.chainMedian` and `densityAnalysis.medians` data member.

        :param dict radiiUpdate: Radii for all atom types.
        :param dict slopesUpdate: Slopes of chain median deviation fraction vs. log b-factor for all atom types.
        :param densityObj: Optional :class:`pdb_eda.ccp4` object.
        :param biopdbObj: Optional `bio.PDB` object.
        :param atomL: Whether or not to calculate statistics for all atoms and assign to `densityAnalysis.atomList`, default as False.
        :param residueL: Whether or not to calculate statistics for all residues and assign to `densityAnalysis.residueList`, default as False.
        :param chainL: Whether or not to calculate statistics for all chains and assign to `densityAnalysis.chainList`, default as False.
        :param recalculate: Whether or not to recalculate if `densityAnalysis.statistics` already exist.

        :type atomL: :py:obj:`True` or :py:obj:`False`
        :type residueL: :py:obj:`True` or :py:obj:`False`
        :type chainL: :py:obj:`True` or :py:obj:`False`
        :type recalculate: :py:obj:`True` or :py:obj:`False`

        :return: :py:obj:`None`
        """
        if self.chainMedian and not recalculate:
            return None
        if not densityObj:
            densityObj = self.densityObj
        if not biopdbObj:
            biopdbObj = self.biopdbObj

        chainClouds = []
        chainPool = []
        chainAvgDensity = []
        residueDict = {}
        chainList = []
        residueList = []
        atomList = []

        currentRadii = {**radiiDefault, **radiiUpdate}
        currentSlopes = {**slopesDefault, **slopesUpdate}
        for residue in biopdbObj.get_residues():
            if residue.id[0] != ' ':
                continue

            residuePool = []
            for atom in residue.child_list:
                resAtom = atom.parent.resname + '_' + atom.name
                if resAtom not in atomType.keys() or atom.get_occupancy() == 0:
                    continue

                ## Calculate atom clouds
                atomClouds = densityObj.findAberrantBlobs(atom.coord, currentRadii[atomType[resAtom]], densityObj.densityCutoff)
                if len(atomClouds) == 0:
                    continue
                elif len(atomClouds) == 1:
                    bestAtomCloud = atomClouds[0]
                else:
                    diffs = [np.linalg.norm(atom.coord - i.centroid) for i in atomClouds]
                    index = diffs.index(min(diffs))
                    bestAtomCloud = atomClouds[index]

                for aCloud in atomClouds:
                    aCloud.atoms = [atom]
                residuePool = residuePool + atomClouds ## For aggregating atom clouds into residue clouds

                atomList.append([residue.parent.id, residue.id[1], atom.parent.resname, atom.name, atomType[resAtom], bestAtomCloud.totalDensity / electrons[resAtom] / atom.get_occupancy(), len(bestAtomCloud.crsList), electrons[resAtom], atom.get_bfactor(), np.linalg.norm(atom.coord - bestAtomCloud.centroid)])
            ## End atom loop

            ## Group connected residue density clouds together from individual atom clouds
            overlap = np.zeros((len(residuePool), len(residuePool)))
            for i in range(len(residuePool)):
                #for j in range(len(residuePool)):
                #    if j <= i:
                #        continue
                for j in range(i+1, len(residuePool)):
                    overlap[i][j] = overlap[j][i] = residuePool[i].testOverlap(residuePool[j])

            resClouds = []
            usedIdx = []
            for i in range(len(residuePool)):
                if i in usedIdx:
                    continue

                newCluster = [n for n, d in enumerate(overlap[i]) if d]
                currCluster = [i] + newCluster
                while len(newCluster):
                    newCluster = {n for x in newCluster for n, d in enumerate(overlap[x]) if n not in currCluster and d}
                    currCluster = currCluster + list(newCluster)

                usedIdx = usedIdx + currCluster
                for idx in currCluster:
                    residuePool[i].merge(residuePool[idx])
                resClouds.append(residuePool[i])

            for cloud in resClouds:
                if len(cloud.atoms) >= 4:
                    if residue.resname in residueDict.keys():
                        residueDict[residue.resname].append(cloud)
                    else:
                        residueDict[residue.resname] = [cloud]

                    totalElectron = sum([electrons[atom.parent.resname + '_' + atom.name] * atom.get_occupancy() for atom in cloud.atoms])
                    residueList.append([residue.parent.id, residue.id[1], residue.resname, cloud.totalDensity / totalElectron, len(cloud.crsList), totalElectron])

            chainPool = chainPool + resClouds ## For aggregating residue clouds into chain clouds
        ## End residue

        ## Group connected chain density clouds together from individual residue clouds
        overlap = np.zeros((len(chainPool), len(chainPool)))
        for i in range(len(chainPool)):
            #for j in range(len(chainPool)):
            #    if j <= i:
            #        continue
            for j in range(i+1, len(chainPool)):
                overlap[i][j] = overlap[j][i] = chainPool[i].testOverlap(chainPool[j])

        usedIdx = []
        for i in range(len(chainPool)):
            if i in usedIdx:
                continue

            newCluster = [n for n, d in enumerate(overlap[i]) if d]
            currCluster = [i] + newCluster
            while len(newCluster):
                newCluster = {n for x in newCluster for n, d in enumerate(overlap[x]) if n not in currCluster and d}
                currCluster = currCluster + list(newCluster)

            usedIdx = usedIdx + currCluster
            for idx in currCluster:
                chainPool[i].merge(chainPool[idx])
            chainClouds.append(chainPool[i])

        for cloud in chainClouds:
            if len(cloud.atoms) <= 50:
                continue
            atom = cloud.atoms[0]
            totalElectron = sum([electrons[atom.parent.resname + '_' + atom.name] * atom.get_occupancy() for atom in cloud.atoms])
            chainList.append([atom.parent.parent.id, atom.parent.id[1], atom.parent.resname, cloud.totalDensity / totalElectron, len(cloud.crsList), totalElectron])
            chainAvgDensity.append(cloud.totalDensity / totalElectron)

        if len(chainAvgDensity) == 0:
            if len(residueList) == 0:
                if len(atomList) == 0:
                    return 0
                else:
                    atomList.sort(key=lambda x: x[5])
                    nVoxel = atomList[len(atomList) // 2][6]
                    totalE = atomList[len(atomList) // 2][7]
                    chainMedian = np.median([x[5] for x in atomList])
            else:
                residueList.sort(key=lambda x: x[3])
                nVoxel = residueList[len(residueList) // 2][4]
                totalE = residueList[len(residueList) // 2][5]
                chainMedian = np.median([x[3] for x in residueList])
        else:
            chainList.sort(key=lambda x: x[3])
            nVoxel = chainList[len(chainList)//2][4]
            totalE = chainList[len(chainList)//2][5]
            chainMedian = np.median(chainAvgDensity)

        # normalize the density by median volume of given atom type
        def normVolumn(row):
            return float(row['density']) / float(row['volume']) * float(medians['volume'][row['atomType']])

        def calcSlope(data):
            ## Less than three data points or all b factors are the same
            if data['chain'].count() <= 2 or all(x == data.iloc[0]['bfactor'] for x in data['bfactor']): 
                return currentSlopes[data.iloc[0]['atomType']]

            slope, intercept, r_vanue, p_value, std_err = stats.linregress(np.log(data['bfactor']), (data['adjDensity']-chainMedian)/chainMedian)
            if p_value > 0.05:
                return currentSlopes[data.iloc[0]['atomType']]
            else:
                return slope

        def getSlope(data):
            return currentSlopes[data.iloc[0]['atomType']]

        def correctFraction(row, slopes, medianBfactor, chainMedian):
            return ((row['adjDensity'] - chainMedian) / chainMedian - (np.log(row['bfactor']) - np.log(medianBfactor.loc[
                medianBfactor.index == row['atomType']])).values * slopes.loc[slopes.index == row['atomType']].values)[0,0]

        try:
            atoms = pandas.DataFrame(atomList, columns=['chain', 'resNum', 'resName', 'atomName', 'atomType', 'density', 'volume', 'electrons', 'bfactor', 'centroidDist'])
            centroidCutoff = atoms['centroidDist'].median() + atoms['centroidDist'].std() * 2
            atoms = atoms[atoms['centroidDist'] < centroidCutoff]  # leave out the atoms that the centroid and atom coordinates are too far away
            medians = atoms.groupby(['atomType']).median()

            ## Normalize by volume
            atoms['adjDensity'] = atoms.apply(lambda row: normVolumn(row), axis=1)
            medians = atoms.groupby(['atomType']).median()
            atoms.loc[atoms.bfactor <= 0, 'bfactor'] = np.nan
            atoms['bfactor'] = atoms.groupby('atomType')['bfactor'].transform(lambda x: x.fillna(x.median()))

            slopes = atoms.groupby('atomType').apply(calcSlope)
            medianBfactor = atoms.groupby('atomType')[['bfactor']].median()

            ## Correct by b-factor
            atoms['chainFraction'] = (atoms['adjDensity'] - chainMedian) / chainMedian
            atoms['correctedFraction'] = atoms.apply(lambda row: correctFraction(row, slopes, medianBfactor, chainMedian), axis=1)
            atoms['correctedDensity'] = atoms['correctedFraction'] * chainMedian + chainMedian
            medians = atoms.groupby(['atomType']).median()
            medians['slopes'] = slopes
        except:
            return 0

        self.chainMedian = chainMedian
        self.chainNvoxel = nVoxel
        self.chainTotalE = totalE
        self.medians = medians
        if atomL:
            self.atomList = atoms
        if residueL:
            self.residueList = residueList
        if chainL:
            self.chainList = chainList


    def getBlobList(self, diffDensityObj=None, recalculate=False):
        """
        Aggregate the difference density map into positive (green) and negative (red) blobs,
        and assign to `densityAnalysis.redBlobList` and `densityAnalysis.greenBlobList`

        :param diffDensityObj: Optional :class:`pdb_eda.ccp4` object.
        :param recalculate: Whether or not to recalculate if `densityAnalysis.statistics` already exist.
        :type recalculate: :py:obj:`True` or :py:obj:`False`

        :return: :py:obj:`None`
        """
        if self.greenBlobList and self.redBlobList and not recalculate:
            return None
        if not diffDensityObj:
            diffDensityObj = self.diffDensityObj

        # find all red/green blobs
        sigma3 = diffDensityObj.diffDensityCutoff

        ## only explore the non-repeating part (<= # xyz intervals) of the density map for blobs
        ncrs = diffDensityObj.header.uniqueNcrs

        ## crs points that are outside 3 sigma
        greenCrsList = np.asarray([[i, j, k] for i in range(ncrs[0]) for j in range(ncrs[1]) for k in range(ncrs[2]) if diffDensityObj.getPointDensityFromCrs([i, j, k]) >= sigma3 ])
        redCrsList = np.asarray([[i, j, k] for i in range(ncrs[0]) for j in range(ncrs[1]) for k in range(ncrs[2]) if diffDensityObj.getPointDensityFromCrs([i, j, k]) <= -sigma3 ])

        ## pairwise distances between all green/red points
        greenDists = scipy.spatial.distance.cdist(greenCrsList, greenCrsList)
        redDists = scipy.spatial.distance.cdist(redCrsList, redCrsList)

        ## group all connected points together into green/red blobs
        dcutoff = np.sqrt(3)  ## the points are considered to be adjacent if one is in the one layer outer box with the other one in the center
        greenBlobList = []
        usedIdx = []
        for i in range(len(greenCrsList)):
            if i in usedIdx:
                continue

            currCluster = [i]
            newCluster = [n for n, d in enumerate(greenDists[i]) if n not in currCluster and d <= dcutoff]
            currCluster = currCluster + newCluster
            while len(newCluster):
                newCluster = {n for x in newCluster for n, d in enumerate(greenDists[x]) if n not in currCluster and d <= dcutoff}
                currCluster = currCluster + list(newCluster)

            usedIdx = usedIdx + currCluster
            blob = ccp4.DensityBlob.fromCrsList([greenCrsList[x] for x in currCluster], diffDensityObj.header, diffDensityObj.density)
            greenBlobList.append(blob)

        redBlobList = []
        usedIdx = []
        for i in range(len(redCrsList)):
            if i in usedIdx:
                continue

            currCluster = [i]
            newCluster = [n for n, d in enumerate(redDists[i]) if n not in currCluster and d <= dcutoff]
            currCluster = currCluster + newCluster
            while len(newCluster):
                newCluster = {n for x in newCluster for n, d in enumerate(redDists[x]) if n not in currCluster and d <= dcutoff}
                currCluster = currCluster + list(newCluster)

            usedIdx = usedIdx + currCluster
            blob = ccp4.DensityBlob.fromCrsList([redCrsList[x] for x in currCluster], diffDensityObj.header, diffDensityObj.density)
            redBlobList.append(blob)

        self.greenBlobList = greenBlobList
        self.redBlobList = redBlobList

    def calcSymmetryAtoms(self, densityObj=None, biopdbObj=None, pdbObj=None, recalculate=False):
        """
        Calculate all the symmetry and nearby cells and keep those have at least on atom within 5 grid points of the non-repeating crs boundary.
        Ref: Biomolecular Crystallography: Principles, Practice, and Application to Structural Biology by Bernhard Rupp.
        Orthogonalization matrix O and deororthogonalization matrix O' are from :class:`pdb_eda.ccp4` object.
        Rotation matrix R and Translation matrix T is from :class:`pdb_eda.pdbParser` object.
        The neighbouring cells can be calculated using formula,
        X' = O(O'(RX + T) + T') = OO'(RX+T) + OT' = RX + T + O[-1/0/1,-1/0/1,-1/0/1].
        Assign the list of :class:`pdb_eda.densityAnalysis.symAtom` instances to `densityAnalysis.symmetryAtoms` data member

        :param densityObj: Optional :class:`pdb_eda.ccp4` object.
        :param biopdbObj: Optional `bio.PDB` object.
        :param pdbObj: Optional :class:`pdb_eda.pdbParser.PDBentry` object.
        :param recalculate: Whether or not to recalculate if `densityAnalysis.statistics` already exist.
        :type recalculate: :py:obj:`True` or :py:obj:`False`

        :return: :py:obj:`None`
        """
        if self.symmetryAtoms and not recalculate:
            return None

        if not densityObj:
            densityObj = self.densityObj
        if not biopdbObj:
            biopdbObj = self.biopdbObj
        if not pdbObj:
            pdbObj = self.pdbObj

        if not self.symmetryAtoms or recalculate:
            pass

        ## For inRangeAtoms, the min/max range of xyz axes (the circumscribed box)
        ncrs = densityObj.header.ncrs
        orginalDensityBox = [densityObj.header.crs2xyzCoord(i) for i in [[c, r, s] for c in [0, ncrs[0]-1] for r in [0, ncrs[1]-1] for s in [0, ncrs[2]-1]]]
        xs = sorted([i[0] for i in orginalDensityBox])
        ys = sorted([i[1] for i in orginalDensityBox])
        zs = sorted([i[2] for i in orginalDensityBox])

        allAtoms = []
        for i in [-1, 0, 1]:
            for j in [-1, 0, 1]:
                for k in [-1, 0, 1]:
                    for r in range(len(pdbObj.header.rotationMats)):
                        if i == 0 and j == 0 and k == 0 and r == 0:
                            inRangeAtoms = list(biopdbObj.get_atoms())
                        else:
                            rMat = pdbObj.header.rotationMats[r]
                            otMat = np.dot(densityObj.header.orthoMat, [i, j, k])
                            atoms = [symAtom(atom) for atom in biopdbObj.get_atoms()]
                            for x in atoms:
                                x.coord = np.dot(rMat[:, 0:3], x.coord) + rMat[:, 3] + otMat

                            ## test if the symmetry atoms are within the range of the original
                            ## convert atom xyz coordinates back to the crs space and check if they are within the original crs range
                            #inRangeAtoms = [x for x in atoms if all([-5 <= densityObj.header.xyz2crsCoord(x.coord)[t] < densityObj.header.uniqueNcrs[t] + 5 for t in range(3)])]

                            inRangeAtoms = [x for x in atoms if xs[0] - 5 <= x.coord[0] <= xs[-1] + 5 and ys[0] - 5 <= x.coord[1] <= ys[-1] + 5 and zs[0] - 5 <= x.coord[2] <= zs[-1] + 5]

                        if len(inRangeAtoms):
                            for x in inRangeAtoms:
                                x.symmetry = [i, j, k, r]
                            allAtoms.extend(inRangeAtoms)

            self.symmetryAtoms = allAtoms

    def calcAtomBlobDists(self, radii={}, slopes={}):
        """
        Calculate `densityAnalysis.symmetryAtoms`, `densityAnalysis.greenBlobList`, `densityAnalysis.redBlobList`, and `densityAnalysis.chainMedian`
        if not already exist, and calculate statistics for positive (green) and negative (red) difference density blobs.

        :return diffMapStats: Difference density map statistics.
        :rtype: :py:obj:`dict`
        """
        if not self.symmetryAtoms:
            self.calcSymmetryAtoms()
        symmetryAtoms = self.symmetryAtoms

        if not self.greenBlobList or not self.redBlobList:
            self.getBlobList()
        greenBlobList = self.greenBlobList
        redBlobList = self.redBlobList

        if not self.chainMedian:
            self.aggregateCloud(radii, slopes)
        chainMedian = self.chainMedian

        ## find the closest atoms to the red/green blobs
        diffMapStats = []
        atomCoords = np.asarray([x.coord for x in symmetryAtoms])
        for blob in greenBlobList + redBlobList:
            ## distanct to the closest atoms
            centroid = np.array(blob.centroid).reshape((1, 3))
            dists = scipy.spatial.distance.cdist(centroid, atomCoords)

            ind = np.argmin(dists[0])
            atom = list(symmetryAtoms)[ind]
            if blob.totalDensity >= 0:
                sign = '+'
            else: sign = '-'
            diffMapStats.append([dists.min(), sign, abs(blob.totalDensity / chainMedian), len(blob.crsList), blob.volume, atom.parent.parent.id, atom.parent.id[1], atom.parent.resname, atom.name, atom.symmetry, atom.coord, blob.centroid])

        return diffMapStats


    def estimateF000(self, densityObj=None, biopdbObj=None, pdbObj=None, recalculate=False):
        """
        Estimate the F000 term as sum of all electrons over the unit cell volume

        :param densityObj: Optional :class:`pdb_eda.ccp4` object.
        :param biopdbObj: Optional `bio.PDB` object.
        :param pdbObj: Optional :class:`pdb_eda.pdbParser.PDBentry` object.
        :param recalculate: Whether or not to recalculate if `densityAnalysis.statistics` already exist.
        :type recalculate: :py:obj:`True` or :py:obj:`False`

        :return: :py:obj:`None`
        """

        if self.f000 and not recalculate:
            return None

        if not densityObj:
            densityObj = self.densityObj
        if not biopdbObj:
            biopdbObj = self.biopdbObj
        if not pdbObj:
            pdbObj = self.pdbObj

        if not self.f000 or recalculate:
            pass

        totalElectrons = 0
        for residue in list(biopdbObj.get_residues()):
            if residue.resname in aaElectrons.keys():
                totalElectrons += aaElectrons[residue.resname]
            else:
                for atom in list(residue.get_atoms()):
                    if atom.name in elementElectrons.keys():
                        totalElectrons += elementElectrons[atom.name]
                totalElectrons += len(list(residue.get_atoms()))  # Add an estimate number of H

        totalElectrons *= len(pdbObj.header.rotationMats)
        asuVolume = densityObj.header.unitVolume * densityObj.header.nintervalX * densityObj.header.nintervalY * densityObj.header.nintervalZ

        self.f000 = totalElectrons/asuVolume


class symAtom:
    """A wrapper class to the `BioPDB.atom` class,
    delegating all BioPDB atom class methods and data members except having its own symmetry and coordination. """

    def __init__(self, atom):
        """
        `pdb_eda.densityAnalysis.symAtom` initializer.

        :param atom: `BioPDB.atom` object.
        """
        self.atom = atom
        self.coord = atom.coord
        self.symmetry = []

    def __getattr__(self, attr):
        return getattr(self.atom, attr)

