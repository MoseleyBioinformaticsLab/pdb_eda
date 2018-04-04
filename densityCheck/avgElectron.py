# !/usr/bin/python3

import ccp4
import pdb as myPDB
#import validationStats

import pandas
import numpy as np
import sys
import os.path
import Bio.PDB as pdb
import matplotlib.pyplot as plt
import scipy.spatial
import copy
import datetime  #; print(str(datetime.datetime.now()))
#from scipy.sparse import coo_matrix
#import crystalContacts

n = 1
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


pdbidfile = sys.argv[1]
pdbids = []
with open(pdbidfile, "r") as fileHandleIn:
    for pdbid in fileHandleIn:
        pdbids.append(pdbid.split(" ; ")[0])

#suffix = sys.argv[2]
#fileHandle = open(sys.argv[2], 'w')
#radii[sys.argv[3]] = float(sys.argv[4])  # for radii optimization
#fileHandleB = open(sys.argv[3], 'w') #for b factor print out

diff = []
for pdbid in pdbids:
    print("working on " + pdbid + ', ', str(datetime.datetime.now()))
    try:
        #atomTypeCount = dict.fromkeys(atomType, 0) ## for atom type composition calculation
        pdbid = pdbid.lower()
        densityObj = ccp4.readFromPDBID(pdbid)
        densityCutoff = np.mean(densityObj.densityArray) + 1.5 * np.std(densityObj.densityArray)

        pdbfile = './pdb/pdb' + pdbid + '.ent'
        if os.path.isfile(pdbfile):
            pass
        else:
            pdbl = pdb.PDBList()
            pdbl.retrieve_pdb_file(pdbid, pdir='./pdb/', file_format="pdb")

        # Bio Python PDB parser
        parser = pdb.PDBParser(QUIET=True)
        structure = parser.get_structure(pdbid, pdbfile)

        ## my own PDB parser
        pdbObj = myPDB.readPDBfile(pdbfile)
        #program = pdbObj.header.program
        #spaceGroup = pdbObj.header.spaceGroup
    except:
        continue

    try:
        diffDensityObj = ccp4.readFromPDBID(pdbid + '_diff')
    except:
        continue

    '''
    ## some RSCC calculation
    valid = validationStats.validationStats(pdbid)
    fo = copy.deepcopy(densityObj)
    fc = copy.deepcopy(densityObj)
    fc.density = densityObj.density - diffDensityObj.density * 2
    #sigma3 = np.mean(densityObj.densityArray) + 1.5 * np.std(densityObj.densityArray)
    sigma3 = 0

    rsccList = valid.getStats(structure, fc, fo, sigma3)
    fileHandle = open("results/rscc." + pdbid + "." + suffix +".txt", 'w')
    print(*rsccList, sep="\n", file=fileHandle)
    fileHandle.close()

    if pdbid == pdbids[len(pdbids)-1]:
        sys.exit()
    else:
        continue
    '''

    ########################################
    ## Aggregate by atom, residue, and chain
    ########################################
    chainClouds = []
    chainAvgDensity = []
    resDict = {}
    chainList = []
    resList = []
    resAvgDensity = []
    atomList = []
    atomAvgDensity = []
    for residue in structure.get_residues():
        if residue.id[0] != ' ': continue

        resClouds = []
        for atom in residue.child_list:
            # if str(atom.parent.id[1]) + '_' + atom.name in contactAtoms: continue

            resAtom = atom.parent.resname + '_' + atom.name
            if resAtom not in atomType.keys():
                continue

            ## Calculate atom blobs
            blobs = densityObj.findAberrantBlobs(atom.coord, radii[atomType[resAtom]], densityCutoff)
            bestBlob = 0
            if len(blobs) == 0:
                continue
            elif len(blobs) == 1:
                bestBlob = blobs[0]
            else:
                diffs = [np.linalg.norm(atom.coord - i.centroid) for i in blobs]
                index = diffs.index(max(diffs))
                bestBlob = blobs[index]

            atomList.append([residue.id[1], resAtom, atomType[resAtom], bestBlob.totalDensity / electrons[resAtom], len(bestBlob.crsList), atom.get_bfactor(), np.linalg.norm(atom.coord - bestBlob.centroid)])
            atomAvgDensity.append(bestBlob.totalDensity / electrons[resAtom])
            #atomTypeCount[resAtom] += 1

            ## Aggregate residue blobs
            for blob in blobs:
                for cloud in resClouds:
                    if cloud.testOverlap(blob):
                        atoms = cloud.atoms
                        cloud.merge(blob)
                        cloud.atoms = atoms + [atom]
                        break
                else:
                    blob.atoms = [atom]
                    resClouds.append(blob)

        for cloud in resClouds:
            if len(cloud.atoms) >= 4:
                if residue.resname in resDict.keys():
                    resDict[residue.resname].append(cloud)
                else:
                    resDict[residue.resname] = [cloud]

                resList.append(str(residue.id[1]) + ', ' + residue.resname + ', ' + str(cloud.totalDensity / sum([electrons[atom.parent.resname + '_' + atom.name] for atom in cloud.atoms])) + ', ' + str(len(cloud.crsList)))
                resAvgDensity.append(cloud.totalDensity / sum([electrons[atom.parent.resname + '_' + atom.name] for atom in cloud.atoms]))

        ## Aggregate chain blobs
        for blob in resClouds:
            for cloud in chainClouds:
                if cloud.testOverlap(blob):
                    atoms = cloud.atoms
                    cloud.merge(blob)
                    cloud.atoms = atoms + blob.atoms
                    break
            else:
                chainClouds.append(blob)

    for cloud in chainClouds:
        if len(cloud.atoms) <= 50: continue
        totalElectron = sum([electrons[atom.parent.resname + '_' + atom.name] for atom in cloud.atoms])
        chainAvgDensity.append(cloud.totalDensity / totalElectron)
        atom = cloud.atoms[0]
        chainList.append(atom.parent.parent.id + ', ' + str(atom.parent.id[1]) + ', ' + atom.parent.resname + ', ' + str(cloud.totalDensity / totalElectron) + ', ' + str(len(cloud.crsList)))

    chainMedian = np.median(chainAvgDensity)

    # reduce atom type to element and single/double/intermediate
    def formatAtomtype(x):
        aa = x.split('_')
        return aa[0] + '_' + aa[1]

    # normalize the density by median volumn of given atom type
    def normVolumn(row):
        return float(row['density']) / float(row['volumn']) * float(medians['volumn'][row['atomType']])

    try:
        atoms = pandas.DataFrame(atomList, columns=['resID', 'resAtom', 'atomType', 'density', 'volumn', 'bfactor', 'centroidDist'])
        centroidCutoff = atoms['centroidDist'].median() + atoms['centroidDist'].std() * 2
        atoms = atoms[atoms['centroidDist'] < centroidCutoff]  # leave out the atoms that the centroid and atom coordinates are too far away

        #atoms['atomType'] = atoms['atomType'].apply(formatAtomtype)
        medians = atoms.groupby(['atomType']).median()
        atoms['adjDensity'] = atoms.apply(lambda row: normVolumn(row), axis=1)
        medians = atoms.groupby(['atomType']).median()
        #medianAdjDen = medians.sort_index()['adjDensity']
    except:
        continue

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

    bfactorMedian = np.median(atoms.bfactor)

    #if n == 1:
      #n = 0
      #print("pdbid", "chainMedian", *sorted(radii.keys()), sep=', ', file=fileHandle)
      #print("pdbid", "chainMedian", *sorted(radii.keys()), sep=', ', file=fileHandleB) ## for print out b factors
      #print("pdbid", *sorted(atomTypeCount.keys()), sep=', ', file=fileHandle) ## for radii optimization, old

    #print(pdbid, chainMedian, *medianAdjDen, sep=", ", file=fileHandle) ## for checking the medians of chain and all atom types
    #print(pdbid, bfactorMedian, *bfactors, sep=", ", file=fileHandleB) ## for print out b factors
    #print(pdbid, *[atomTypeCount[key] for key in sorted(atomTypeCount.keys())], sep=', ', file=fileHandle) ## for atom type composition

    #diff.append((chainMedian - medianAdjDen[int(sys.argv[5])]) / chainMedian)  ## for radii optimization


    ###########################
    # find all red/green blobs
    ###########################
    sigma3 = np.mean(diffDensityObj.densityArray) + 3 * np.std(diffDensityObj.densityArray)

    ## only explore the non-repeating part (<= # xyz intervals) of the density map for blobs
    ncrs = [i for i in diffDensityObj.header.ncrs]
    if diffDensityObj.header.xyzInterval[diffDensityObj.header.col2xyz - 1] < ncrs[0]:
        ncrs[0] = diffDensityObj.header.xyzInterval[diffDensityObj.header.col2xyz - 1]
    if diffDensityObj.header.xyzInterval[diffDensityObj.header.row2xyz - 1] < ncrs[1]:
        ncrs[1] = diffDensityObj.header.xyzInterval[diffDensityObj.header.row2xyz - 1]
    if diffDensityObj.header.xyzInterval[diffDensityObj.header.sec2xyz - 1] < ncrs[2]:
        ncrs[2] = diffDensityObj.header.xyzInterval[diffDensityObj.header.sec2xyz - 1]

    ## crs points that are outside 3 sigma
    greenCrsList = np.asarray([[i, j, k] for i in range(ncrs[0]) for j in range(ncrs[1]) for k in range(ncrs[2]) if diffDensityObj.getPointDensityFromCrs([i, j, k]) >= sigma3 > 0])
    redCrsList = np.asarray([[i, j, k] for i in range(ncrs[0]) for j in range(ncrs[1]) for k in range(ncrs[2]) if diffDensityObj.getPointDensityFromCrs([i, j, k]) <= -sigma3 < 0])

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


    '''
    ## For inRangeAtoms0, the min/max range of xyz axes (the circumscribed box)
    orginalDensityBox = [diffDensityObj.header.crs2xyzCoord(i) for i in [[c, r, s] for c in [0, ncrs[0]-1] for r in [0, ncrs[1]-1] for s in [0, ncrs[2]-1]]]
    xs = sorted([i[0] for i in orginalDensityBox])
    ys = sorted([i[1] for i in orginalDensityBox])
    zs = sorted([i[2] for i in orginalDensityBox])

    ## For inRangeAtoms00, the mathematical calculation of checking if a point is inside a parallelepiped box
    ## https://math.stackexchange.com/questions/1472049/check-if-a-point-is-inside-a-rectangular-shaped-area-3d
    p1 = diffDensityObj.header.crs2xyzCoord([0, 0, 0])
    p2 = diffDensityObj.header.crs2xyzCoord([ncrs[0]-1, 0, 0])
    p4 = diffDensityObj.header.crs2xyzCoord([0, ncrs[1]-1, 0])
    p5 = diffDensityObj.header.crs2xyzCoord([0, 0, ncrs[2]-1])
    u = np.cross(p1 - p4, p1 - p5)
    v = np.cross(p1 - p2, p1 - p5)
    w = np.cross(p1 - p2, p1 - p4)
    up1 = round(np.dot(u, p1), 6)
    up2 = round(np.dot(u, p2), 6)
    vp1 = round(np.dot(v, p1), 6)
    vp4 = round(np.dot(v, p4), 6)
    wp1 = round(np.dot(w, p1), 6)
    wp5 = round(np.dot(w, p5), 6)
    '''

    ## calculate all symmetry and nearby cells
    ## Biomolecular Crystallography: Principles, Practice, and Application to Structural Biology by Bernhard Rupp
    ## Orthogonalization matrix O and deororthogonalization matrix O' are from 'ccp4'
    ## Rotation Matrix is from 'myPDB'
    ## The neighbering cells can be calculated using formula, X' = O(O'(RX + T) + T') = OO'(RX+T) + OT' = RX + T + O[-1/0/1,-1/0/1,-1/0/1]
    atomCoords = np.array([i.coord for i in structure.get_atoms()])
    allAtoms = []
    allAtoms.extend(atomCoords)

    #print(str(datetime.now()))
    allAtoms = []
    for i in [-1, 0, 1]:
        for j in [-1, 0, 1]:
            for k in [-1, 0, 1]:
                for r in range(len(pdbObj.header.rotationMats)):
                    if i == 0 and j == 0 and k == 0 and r == 0:
                        inRangeAtoms = list(structure.get_atoms())
                    else:
                        rMat = pdbObj.header.rotationMats[r]
                        otMat = np.dot(densityObj.header.orthoMat, [i, j, k])
                        atoms = copy.deepcopy(list(structure.get_atoms()))
                        for x in atoms:
                            x.coord = np.dot(rMat[:, 0:3], x.coord) + rMat[:, 3] + otMat

                        ## test if the symmetry atoms are within the range of the original
                        ## convert atom xyz coordinates back to the crs space and check if they are within the original crs range
                        inRangeAtoms = [x for x in atoms if all([-5 <= diffDensityObj.header.xyz2crsCoord(x.coord)[t] < ncrs[t] + 5 for t in range(3)])]

                        ## other methods that can somewhat validate the above
                        # inRangeAtomsBox = [x for x in symAtoms if xs[0] - 5 <= x[0] <= xs[-1] + 5 and ys[0] - 5 <= x[1] <= ys[-1] + 5 and zs[0] - 5 <= x[2] <= zs[-1] + 5]
                        # inRangeAtomsMath = [x for x in symAtoms if (up1 >= round(np.dot(u, x),6) >= up2 or up1 <= round(np.dot(u, x),6) <= up2) and (vp1 >= round(np.dot(v, x),6) >= vp4 or vp1 <= round(np.dot(v, x),6) <= vp4) and (wp1 >= round(np.dot(w, x),6) >= wp5 or wp1 <= round(np.dot(w, x),6) <= wp5)]

                    if len(inRangeAtoms):
                        for x in inRangeAtoms:
                            x.symmetry = [i, j, k, r]
                        allAtoms.extend(inRangeAtoms)

    atomCoords = np.asarray([x.coord for x in allAtoms])

    ## find nearby atoms to the red/green blobs
    diffMapStats = []
    isolatedBlobs = []
    for blob in greenBlobList + redBlobList:
        ## distanct to the closest atoms
        centroid = np.array(blob.centroid).reshape((1, 3))
        dists = scipy.spatial.distance.cdist(centroid, atomCoords)

        ind = np.argmin(dists[0])
        atom = list(allAtoms)[ind]
        if blob.totalDensity >= 0:
            sign = '+'
        else: sign = '-'
        diffMapStats.append([dists.min(), sign, abs(blob.totalDensity / chainMedian), blob.volume, atom.parent.parent.id, atom.parent.id[1], atom.parent.resname, atom.name, atom.symmetry, atom.coord, blob.centroid])


    dists = [i[0] for i in diffMapStats]
    histogram = plt.hist(dists, 200)
    plt.savefig('../atom-red-avg-distance/' + pdbid + '.png')
    plt.close()

    #diffMapStats.sort(key=lambda x: x[2], reverse=True)  # sort by number of electron
    #diffMapStats.sort(key=lambda x: x[4])  # sort by distance

    #dists = [row[4] for row in diffMapStats]
    #plt.hist(dists, bins=200)

    #dists = [row[6] for row in diffMapStats]
    #plt.hist(dists, bins=200)
    #plt.savefig('../' + pdbid + '.png')
    #plt.close()

    #print(len([i for i in dists if i <= 2]))
    #print(len(dists))

    '''
    dists = [i[1] for i in diffMapStats]
    histogram = plt.hist(dists, bins=np.arange(min(dists), max(dists) + 0.02, 0.02))
    plt.close()
    mode = histogram[1][np.argmax(histogram[0][0:250]) + 1]  # maximum bin < 5A (5/0.02 = 250) distance
    logmode = np.log(mode)

    logdists = [np.log(i) for i in dists]
    leftside = [i for i in logdists if i < logmode]
    dev1 = np.sqrt(sum([(i - logmode) ** 2 for i in leftside]) / len(leftside))

    bothside = [i for i in logdists if i < logmode + 2 * dev1]
    cutoff = np.mean(bothside) + 2 * np.std(bothside)

    plt.hist(logdists, 200)
    plt.axvline(x=logmode, color='red')
    plt.axvline(x=logmode + dev1, color='orange')
    plt.axvline(x=logmode + 2 * dev1, color='yellow')
    plt.axvline(x=cutoff, color='green')
    plt.savefig('../atom-blue-distance/' + pdbid + '.log.png')
    plt.close()

    model = [i for i in dists if i < 10]
    plt.hist(model, 200)
    plt.axvline(x=mode, color='red')
    plt.axvline(x=np.exp(logmode + dev1), color='orange')
    plt.axvline(x=np.exp(logmode + 2 * dev1), color='yellow')
    plt.axvline(x=np.exp(cutoff), color='green')
    plt.savefig('../atom-blue-distance/' + pdbid + '.original.png')
    plt.close()

    model = [i for i in diffMapStats if i[1] < np.exp(cutoff)]

    plt.hist([i[2] for i in model], bins=100)
    plt.axvline(x=np.exp(cutoff), color='green')
    plt.savefig('../min-red-blue-distance/' + pdbid + '.2.png')
    plt.close()

    plt.hist([i[3] for i in model], bins=100)
    plt.axvline(x=np.exp(cutoff), color='green')
    plt.savefig('../red-atom-avg-distance/' + pdbid + '.2.png')
    plt.close()

    plt.hist([i[4] for i in model], bins=100)
    plt.axvline(x=np.exp(cutoff), color='green')
    plt.savefig('../min-red-atom-distance/' + pdbid + '.2.png')
    plt.close()
    '''

#print(np.nanmean(diff), np.nanmedian(diff), file=fileHandle) ## for radii optimization

#fileHandle = open("results/cen." + pdbid + ".txt", 'w') ## form print out centroid-coordinates difference
#print(*atomList, sep='\n', file=fileHandle) ## for single atoms

#fileHandle.close()
#fileHandleB.close()

print('Done! ', str(datetime.datetime.now()))

