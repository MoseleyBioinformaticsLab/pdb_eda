# !/usr/bin/python3
"""
densityAnalysis.py
    Reads in pdbid and has functions to analyze its electron density.
"""

import copy
import urllib.request
import datetime  #print(str(datetime.datetime.now()))
import os.path

import pandas
import numpy as np
import Bio.PDB as biopdb
import scipy.spatial
from scipy import stats

from . import ccp4
from . import pdbParser
from . import validationStats


elementElectron = {'C': 6, 'N': 7, 'O': 8, 'P': 15, 'S': 16}
electrons = {'GLY_N': 8, 'GLY_CA': 8, 'GLY_C': 6, 'GLY_O': 8, 'GLY_OXT': 8,
             'ALA_N': 8, 'ALA_CA': 7, 'ALA_C': 6, 'ALA_O': 8, 'ALA_CB': 9, 'ALA_OXT': 8,
             'VAL_N': 8, 'VAL_CA': 7, 'VAL_C': 6, 'VAL_O': 8, 'VAL_CB': 7, 'VAL_CG1': 9, 'VAL_CG2': 9, 'VAL_OXT': 8,
             'LEU_N': 8, 'LEU_CA': 7, 'LEU_C': 6, 'LEU_O': 8, 'LEU_CB': 8, 'LEU_CG': 7, 'LEU_CD1': 9, 'LEU_CD2': 9, 'LEU_OXT': 8,
             'ILE_N': 8, 'ILE_CA': 7, 'ILE_C': 6, 'ILE_O': 8, 'ILE_CB': 7, 'ILE_CG1': 8, 'ILE_CG2': 9, 'ILE_CD1': 9, 'ILE_OXT': 8,
             'MET_N': 8, 'MET_CA': 7, 'MET_C': 6, 'MET_O': 8, 'MET_CB': 8, 'MET_CG': 8, 'MET_SD': 16, 'MET_CE': 9, 'MET_OXT': 8,
             'PHE_N': 8, 'PHE_CA': 7, 'PHE_C': 6, 'PHE_O': 8, 'PHE_CB': 8, 'PHE_CG': 6, 'PHE_CD1': 7, 'PHE_CD2': 7, 'PHE_CE1': 7, 'PHE_CE2': 7, 'PHE_CZ': 7, 'PHE_OXT': 8,
             'TRP_N': 8, 'TRP_CA': 7, 'TRP_C': 6, 'TRP_O': 8, 'TRP_CB': 8, 'TRP_CG': 6, 'TRP_CD1': 7, 'TRP_CD2': 6, 'TRP_NE1': 8, 'TRP_CE2': 6, 'TRP_CE3': 7, 'TRP_CZ2': 7, 'TRP_CZ3': 7, 'TRP_CH2': 7, 'TRP_OXT': 8,
             'PRO_N': 7, 'PRO_CA': 7, 'PRO_C': 6, 'PRO_O': 8, 'PRO_CB': 8, 'PRO_CG': 8, 'PRO_CD': 8, 'PRO_OXT': 8,
             'SER_N': 8, 'SER_CA': 7, 'SER_C': 6, 'SER_O': 8, 'SER_CB': 8, 'SER_OG': 9, 'SER_OXT': 8,
             'THR_N': 8, 'THR_CA': 7, 'THR_C': 6, 'THR_O': 8, 'THR_CB': 7, 'THR_OG1': 9, 'THR_CG2': 9, 'THR_OXT': 8,
             'CYS_N': 8, 'CYS_CA': 7, 'CYS_C': 6, 'CYS_O': 8, 'CYS_CB': 8, 'CYS_SG': 17, 'CYS_OXT': 8,
             'TYR_N': 8, 'TYR_CA': 7, 'TYR_C': 6, 'TYR_O': 8, 'TYR_CB': 8, 'TYR_CG': 6, 'TYR_CD1': 7, 'TYR_CD2': 7, 'TYR_CE1': 7, 'TYR_CE2': 7, 'TYR_CZ': 6, 'TYR_OH': 9, 'TYR_OXT': 8,
             'ASN_N': 8, 'ASN_CA': 7, 'ASN_C': 6, 'ASN_O': 8, 'ASN_CB': 8, 'ASN_CG': 6, 'ASN_OD1': 8, 'ASN_ND2': 9, 'ASN_OXT': 8,
             'GLN_N': 8, 'GLN_CA': 7, 'GLN_C': 6, 'GLN_O': 8, 'GLN_CB': 8, 'GLN_CG': 8, 'GLN_CD': 6, 'GLN_OE1': 8, 'GLN_NE2': 9, 'GLN_OXT': 8,
             'ASP_N': 8, 'ASP_CA': 7, 'ASP_C': 6, 'ASP_O': 8, 'ASP_CB': 8, 'ASP_CG': 6, 'ASP_OD1': 8, 'ASP_OD2': 8, 'ASP_OXT': 8,
             'GLU_N': 8, 'GLU_CA': 7, 'GLU_C': 6, 'GLU_O': 8, 'GLU_CB': 8, 'GLU_CG': 8, 'GLU_CD': 6, 'GLU_OE1': 8, 'GLU_OE2': 8, 'GLU_OXT': 8,
             'LYS_N': 8, 'LYS_CA': 7, 'LYS_C': 6, 'LYS_O': 8, 'LYS_CB': 8, 'LYS_CG': 8, 'LYS_CD': 8, 'LYS_CE': 8, 'LYS_NZ': 9, 'LYS_OXT': 8,
             'ARG_N': 8, 'ARG_CA': 7, 'ARG_C': 6, 'ARG_O': 8, 'ARG_CB': 8, 'ARG_CG': 8, 'ARG_CD': 8, 'ARG_NE': 8, 'ARG_CZ': 6, 'ARG_NH1': 8, 'ARG_NH2': 8, 'ARG_OXT': 8,
             'HIS_N': 8, 'HIS_CA': 7, 'HIS_C': 6, 'HIS_O': 8, 'HIS_CB': 8, 'HIS_CG': 6, 'HIS_ND1': 7, 'HIS_CD2': 7, 'HIS_CE1': 7, 'HIS_NE2': 8, 'HIS_OXT': 8}

atomType = {'GLY_N': 'N_single_bb', 'GLY_CA': 'C_single_bb', 'GLY_C': 'C_double_bb', 'GLY_O': 'O_double_bb', 'GLY_OXT': 'O_intermediate',
            'ALA_N': 'N_single_bb', 'ALA_CA': 'C_single_bb', 'ALA_C': 'C_double_bb', 'ALA_O': 'O_double_bb', 'ALA_CB': 'C_single', 'ALA_OXT': 'O_intermediate',
            'VAL_N': 'N_single_bb', 'VAL_CA': 'C_single_bb', 'VAL_C': 'C_double_bb', 'VAL_O': 'O_double_bb', 'VAL_CB': 'C_single', 'VAL_CG1': 'C_single', 'VAL_CG2': 'C_single', 'VAL_OXT': 'O_intermediate',
            'LEU_N': 'N_single_bb', 'LEU_CA': 'C_single_bb', 'LEU_C': 'C_double_bb', 'LEU_O': 'O_double_bb', 'LEU_CB': 'C_single', 'LEU_CG': 'C_single', 'LEU_CD1': 'C_single', 'LEU_CD2': 'C_single', 'LEU_OXT': 'O_intermediate',
            'ILE_N': 'N_single_bb', 'ILE_CA': 'C_single_bb', 'ILE_C': 'C_double_bb', 'ILE_O': 'O_double_bb', 'ILE_CB': 'C_single', 'ILE_CG1': 'C_single', 'ILE_CG2': 'C_single', 'ILE_CD1': 'C_single', 'ILE_OXT': 'O_intermediate',
            'MET_N': 'N_single_bb', 'MET_CA': 'C_single_bb', 'MET_C': 'C_double_bb', 'MET_O': 'O_double_bb', 'MET_CB': 'C_single', 'MET_CG': 'C_single', 'MET_SD': 'S_single', 'MET_CE': 'C_single', 'MET_OXT': 'O_intermediate',
            'PHE_N': 'N_single_bb', 'PHE_CA': 'C_single_bb', 'PHE_C': 'C_double_bb', 'PHE_O': 'O_double_bb', 'PHE_CB': 'C_single', 'PHE_CG': 'C_intermediate', 'PHE_CD1': 'C_intermediate', 'PHE_CD2': 'C_intermediate', 'PHE_CE1': 'C_intermediate', 'PHE_CE2': 'C_intermediate', 'PHE_CZ': 'C_intermediate', 'PHE_OXT': 'O_intermediate',
            'TRP_N': 'N_single_bb', 'TRP_CA': 'C_single_bb', 'TRP_C': 'C_double_bb', 'TRP_O': 'O_double_bb', 'TRP_CB': 'C_single', 'TRP_CG': 'C_intermediate', 'TRP_CD1': 'C_intermediate', 'TRP_CD2': 'C_intermediate', 'TRP_NE1': 'N_intermediate', 'TRP_CE2': 'C_intermediate', 'TRP_CE3': 'C_intermediate', 'TRP_CZ2': 'C_intermediate', 'TRP_CZ3': 'C_intermediate', 'TRP_CH2': 'C_intermediate', 'TRP_OXT': 'O_intermediate',
            'PRO_N': 'N_single_bb', 'PRO_CA': 'C_single_bb', 'PRO_C': 'C_double_bb', 'PRO_O': 'O_double_bb', 'PRO_CB': 'C_single', 'PRO_CG': 'C_single', 'PRO_CD': 'C_single', 'PRO_OXT': 'O_intermediate',
            'SER_N': 'N_single_bb', 'SER_CA': 'C_single_bb', 'SER_C': 'C_double_bb', 'SER_O': 'O_double_bb', 'SER_CB': 'C_single', 'SER_OG': 'O_single', 'SER_OXT': 'O_intermediate',
            'THR_N': 'N_single_bb', 'THR_CA': 'C_single_bb', 'THR_C': 'C_double_bb', 'THR_O': 'O_double_bb', 'THR_CB': 'C_single', 'THR_OG1': 'O_single', 'THR_CG2': 'C_single', 'THR_OXT': 'O_intermediate',
            'CYS_N': 'N_single_bb', 'CYS_CA': 'C_single_bb', 'CYS_C': 'C_double_bb', 'CYS_O': 'O_double_bb', 'CYS_CB': 'C_single', 'CYS_SG': 'S_single', 'CYS_OXT': 'O_intermediate',
            'TYR_N': 'N_single_bb', 'TYR_CA': 'C_single_bb', 'TYR_C': 'C_double_bb', 'TYR_O': 'O_double_bb', 'TYR_CB': 'C_single', 'TYR_CG': 'C_intermediate', 'TYR_CD1': 'C_intermediate', 'TYR_CD2': 'C_intermediate', 'TYR_CE1': 'C_intermediate', 'TYR_CE2': 'C_intermediate', 'TYR_CZ': 'C_intermediate', 'TYR_OH': 'O_single', 'TYR_OXT': 'O_intermediate',
            'ASN_N': 'N_single_bb', 'ASN_CA': 'C_single_bb', 'ASN_C': 'C_double_bb', 'ASN_O': 'O_double_bb', 'ASN_CB': 'C_single', 'ASN_CG': 'C_double', 'ASN_OD1': 'O_double', 'ASN_ND2': 'N_single', 'ASN_OXT': 'O_intermediate',
            'GLN_N': 'N_single_bb', 'GLN_CA': 'C_single_bb', 'GLN_C': 'C_double_bb', 'GLN_O': 'O_double_bb', 'GLN_CB': 'C_single', 'GLN_CG': 'C_single', 'GLN_CD': 'C_double', 'GLN_OE1': 'O_double', 'GLN_NE2': 'N_single', 'GLN_OXT': 'O_intermediate',
            'ASP_N': 'N_single_bb', 'ASP_CA': 'C_single_bb', 'ASP_C': 'C_double_bb', 'ASP_O': 'O_double_bb', 'ASP_CB': 'C_single', 'ASP_CG': 'C_double', 'ASP_OD1': 'O_intermediate', 'ASP_OD2': 'O_intermediate', 'ASP_OXT': 'O_intermediate',
            'GLU_N': 'N_single_bb', 'GLU_CA': 'C_single_bb', 'GLU_C': 'C_double_bb', 'GLU_O': 'O_double_bb', 'GLU_CB': 'C_single', 'GLU_CG': 'C_single', 'GLU_CD': 'C_double', 'GLU_OE1': 'O_intermediate', 'GLU_OE2': 'O_intermediate', 'GLU_OXT': 'O_intermediate',
            'LYS_N': 'N_single_bb', 'LYS_CA': 'C_single_bb', 'LYS_C': 'C_double_bb', 'LYS_O': 'O_double_bb', 'LYS_CB': 'C_single', 'LYS_CG': 'C_single', 'LYS_CD': 'C_single', 'LYS_CE': 'C_single', 'LYS_NZ': 'N_single', 'LYS_OXT': 'O_intermediate',
            'ARG_N': 'N_single_bb', 'ARG_CA': 'C_single_bb', 'ARG_C': 'C_double_bb', 'ARG_O': 'O_double_bb', 'ARG_CB': 'C_single', 'ARG_CG': 'C_single', 'ARG_CD': 'C_single', 'ARG_NE': 'N_intermediate', 'ARG_CZ': 'C_double', 'ARG_NH1': 'N_intermediate', 'ARG_NH2': 'N_intermediate', 'ARG_OXT': 'O_intermediate',
            'HIS_N': 'N_single_bb', 'HIS_CA': 'C_single_bb', 'HIS_C': 'C_double_bb', 'HIS_O': 'O_double_bb', 'HIS_CB': 'C_single', 'HIS_CG': 'C_intermediate', 'HIS_ND1': 'N_intermediate', 'HIS_CD2': 'C_intermediate', 'HIS_CE1': 'C_intermediate', 'HIS_NE2': 'N_intermediate', 'HIS_OXT': 'O_intermediate'}

## Data originally from https://arxiv.org/pdf/0804.2488.pdf
radii = {'C_single': 0.91, 'C_double': 0.71, 'C_intermediate': 0.76, 'C_single_bb': 0.74, 'C_double_bb': 0.64,
         'O_single': 0.88, 'O_double': 0.83, 'O_intermediate': 0.91, 'O_double_bb': 0.75,
         'N_single': 1.03, 'N_intermediate': 0.83, 'N_single_bb': 0.73, #'N_double': 0.74,
         'S_single': 0.80}

ccp4urlPrefix = "http://www.ebi.ac.uk/pdbe/coordinates/files/"
ccp4urlSuffix = ".ccp4"
ccp4folder = './ccp4_data/'
pdbfolder = './pdb_data/'


def fromPDBid(pdbid, ccp4density=True, ccp4diff=True, pdbbio=True, pdbi=True, downloadFile=True):
    """RETURNS DensityAnalysis object given the PARAMETER pdbid."""
    pdbid = pdbid.lower()
    print("working on " + pdbid + ', ', str(datetime.datetime.now()))

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
    def __init__(self, pdbid, densityObj=None, diffDensityObj=None, biopdbObj=None, pdbObj=None):
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


    def validation(self, densityObj=None, diffDensityObj=None, biopdbObj=None, recalculate=False):
        """
        RETURNS
            validation statistics (RSR and RSCC) given,
        PARAMS
            densityObj: density object, use initialized data member if not provided
            diffDensityObj: difference density object, use initialized data member if not provided
            biopdbObj: BioPDB object, use initialized data member if not provided
            recalculate: recalculate if already exist, default as False
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


    def aggregateCloud(self, densityObj=None, biopdbObj=None, atomL=False, residueL=False, chainL=False, recalculate=False):
        """
        RETURNS
            chainMedian and medians (by atom type) data member given,
        PARAMS
            densityObj: density object, use initialized data member if not provided
            biopdbObj: BioPDB object, use initialized data member if not provided
            atomL: whether or not return a full atom list, default as False
            residueL: whether or not return a full residue list, default as False
            chainL: whether or not return a full chain list, default as False
            recalculate: recalculate if already exist, default as False
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
        for residue in biopdbObj.get_residues():
            if residue.id[0] != ' ':
                continue

            residuePool = []
            for atom in residue.child_list:
                resAtom = atom.parent.resname + '_' + atom.name
                if resAtom not in atomType.keys():
                    continue

                ## Calculate atom clouds
                atomClouds = densityObj.findAberrantBlobs(atom.coord, radii[atomType[resAtom]], densityObj.densityCutoff)
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

                atomList.append([residue.parent.id, residue.id[1], atom.parent.resname, atom.name, atomType[resAtom], bestAtomCloud.totalDensity / electrons[resAtom], len(bestAtomCloud.crsList), atom.get_bfactor(), np.linalg.norm(atom.coord - bestAtomCloud.centroid)])
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

                    residueList.append([residue.parent.id, residue.id[1], residue.resname, cloud.totalDensity / sum([electrons[atom.parent.resname + '_' + atom.name] for atom in cloud.atoms]), len(cloud.crsList)])

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
            totalElectron = sum([electrons[atom.parent.resname + '_' + atom.name] for atom in cloud.atoms])
            chainList.append([atom.parent.parent.id, atom.parent.id[1], atom.parent.resname, cloud.totalDensity / totalElectron, len(cloud.crsList)])
            chainAvgDensity.append(cloud.totalDensity / totalElectron)

        chainMedian = np.median(chainAvgDensity)

        # normalize the density by median volumn of given atom type
        def normVolumn(row):
            return float(row['density']) / float(row['volumn']) * float(medians['volumn'][row['atomType']])

        try:
            atoms = pandas.DataFrame(atomList, columns=['chain', 'resNum', 'resName', 'atomName', 'atomType', 'density', 'volumn', 'bfactor', 'centroidDist'])
            centroidCutoff = atoms['centroidDist'].median() + atoms['centroidDist'].std() * 2
            atoms = atoms[atoms['centroidDist'] < centroidCutoff]  # leave out the atoms that the centroid and atom coordinates are too far away

            medians = atoms.groupby(['atomType']).median()
            atoms['adjDensity'] = atoms.apply(lambda row: normVolumn(row), axis=1)
            medians = atoms.groupby(['atomType']).median()
            #medianAdjDen = medians.sort_index()['adjDensity']
        except:
            return 0

        medianAdjDen = []
        bfactors = []
        for key in sorted(radii.keys()):
            try:
                medianAdjDen.append(medians['adjDensity'][key])
            except:
                medianAdjDen.append(np.nan)

            ## b factors calculation
            try:
                bfactors.append(medians['bfactor'][key])
            except:
                bfactors.append(np.nan)

        self.chainMedian = chainMedian
        self.medians = medians
        if atomL:
            self.atomList = atomList
        if residueL:
            self.residueList = residueList
        if chainL:
            self.chainList = chainList


    def getBlobList(self, diffDensityObj=None, recalculate=False):
        """
        RETURNS
            green and red blob lists as data member given,
        PARAMS
            diffDensityObj: difference density object, use initialized data member if not provided
            recalculate: recalculate if already exist, default as False
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
        RETURNS
            symmetryAtoms as data member given,
        PARAMS
            densityObj: density object, use initialized data member if not provided
            diffDensityObj: difference density object, use initialized data member if not provided
            pdbObj: my pdb object, use initialized data member if not provided
            recalculate: recalculate if already exist, default as False

        ## calculate all symmetry and nearby cells with in 5 grid points of the non-repeating crs boundary
        ## Biomolecular Crystallography: Principles, Practice, and Application to Structural Biology by Bernhard Rupp
        ## Orthogonalization matrix O and deororthogonalization matrix O' are from 'ccp4'
        ## Rotation Matrix is from 'myPDB'
        ## The neighbering cells can be calculated using formula, X' = O(O'(RX + T) + T') = OO'(RX+T) + OT' = RX + T + O[-1/0/1,-1/0/1,-1/0/1]
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

    def calcAtomBlobDists(self):
        """
        RETURNS
            diffMapStats given,
        PARAMS
            symmetryAtoms, greenBlobList, redBlobList:, chainMedian: calculate not exist
        """
        if not self.symmetryAtoms:
            self.calcSymmetryAtoms()
        symmetryAtoms = self.symmetryAtoms

        if not self.greenBlobList or not self.redBlobList:
            self.getBlobList()
        greenBlobList = self.greenBlobList
        redBlobList = self.redBlobList

        if not self.chainMedian:
            self.aggregateCloud()
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


class symAtom:
    """ A wrapper class to the BioPDB atom class, delegate all BioPDB atom class method and data member except having its own symmetry and coordination """
    def __init__(self, atom):
        self.atom = atom
        self.coord = atom.coord
        self.symmetry = []

    def __getattr__(self, attr):
        return getattr(self.atom, attr)

