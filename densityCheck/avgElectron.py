# !/usr/bin/python3

import pandas
import ccp4
import numpy as np
import sys
import os.path
import pdb as myPDB
import Bio.PDB as pdb
#import crystalContacts

elementElectron = {'C': 6, 'N': 7, 'O': 8, 'P': 15, 'S': 16}
electrons = {'GLY_N': 8, 'GLY_CA': 8, 'GLY_C': 6, 'GLY_O': 8, 'GLY_OXT': 9,
             'ALA_N': 8, 'ALA_CA': 7, 'ALA_C': 6, 'ALA_O': 8, 'ALA_CB': 9, 'ALA_OXT': 9,
             'VAL_N': 8, 'VAL_CA': 7, 'VAL_C': 6, 'VAL_O': 8, 'VAL_CB': 7, 'VAL_CG1': 9, 'VAL_CG2': 9, 'VAL_OXT': 9,
             'LEU_N': 8, 'LEU_CA': 7, 'LEU_C': 6, 'LEU_O': 8, 'LEU_CB': 8, 'LEU_CG': 8, 'LEU_CD1': 9, 'LEU_CD2': 9, 'LEU_OXT': 9,
             'ILE_N': 8, 'ILE_CA': 7, 'ILE_C': 6, 'ILE_O': 8, 'ILE_CB': 7, 'ILE_CG1': 8, 'ILE_CG2': 9, 'ILE_CD1': 9, 'ILE_OXT': 9,
             'MET_N': 8, 'MET_CA': 7, 'MET_C': 6, 'MET_O': 8, 'MET_CB': 8, 'MET_CG': 8, 'MET_SD': 16, 'MET_CE': 9, 'MET_OXT': 9,
             'PHE_N': 8, 'PHE_CA': 7, 'PHE_C': 6, 'PHE_O': 8, 'PHE_CB': 8, 'PHE_CG': 6, 'PHE_CD1': 7, 'PHE_CD2': 7, 'PHE_CE1': 7, 'PHE_CE2': 7, 'PHE_CZ': 7, 'PHE_OXT': 9,
             'TRP_N': 8, 'TRP_CA': 7, 'TRP_C': 6, 'TRP_O': 8, 'TRP_CB': 8, 'TRP_CG': 6, 'TRP_CD1': 7, 'TRP_CD2': 6, 'TRP_NE1': 8, 'TRP_CE2': 6, 'TRP_CE3': 7, 'TRP_CZ2': 7, 'TRP_CZ3': 7, 'TRP_CH2': 7, 'TRP_OXT': 9,
             'PRO_N': 7, 'PRO_CA': 7, 'PRO_C': 6, 'PRO_O': 8, 'PRO_CB': 8, 'PRO_CG': 8, 'PRO_CD': 8, 'PRO_OXT': 9,
             'SER_N': 8, 'SER_CA': 7, 'SER_C': 6, 'SER_O': 8, 'SER_CB': 8, 'SER_OG': 9, 'SER_OXT': 9,
             'THR_N': 8, 'THR_CA': 7, 'THR_C': 6, 'THR_O': 8, 'THR_CB': 7, 'THR_OG1': 9, 'THR_CG2': 9, 'THR_OXT': 9,
             'CYS_N': 8, 'CYS_CA': 7, 'CYS_C': 6, 'CYS_O': 8, 'CYS_CB': 8, 'CYS_SG': 17, 'CYS_OXT': 9,
             'TYR_N': 8, 'TYR_CA': 7, 'TYR_C': 6, 'TYR_O': 8, 'TYR_CB': 8, 'TYR_CG': 6, 'TYR_CD1': 7, 'TYR_CD2': 7, 'TYR_CE1': 7, 'TYR_CE2': 7, 'TYR_CZ': 6, 'TYR_OH': 9, 'TYR_OXT': 9,
             'ASN_N': 8, 'ASN_CA': 7, 'ASN_C': 6, 'ASN_O': 8, 'ASN_CB': 8, 'ASN_CG': 6, 'ASN_OD1': 8, 'ASN_ND2': 9, 'ASN_OXT': 9,
             'GLN_N': 8, 'GLN_CA': 7, 'GLN_C': 6, 'GLN_O': 8, 'GLN_CB': 8, 'GLN_CG': 8, 'GLN_CD': 6, 'GLN_OE1': 8, 'GLN_NE2': 9, 'GLN_OXT': 9,
             'ASP_N': 8, 'ASP_CA': 7, 'ASP_C': 6, 'ASP_O': 8, 'ASP_CB': 8, 'ASP_CG': 6, 'ASP_OD1': 8, 'ASP_OD2': 9, 'ASP_OXT': 9,
             'GLU_N': 8, 'GLU_CA': 7, 'GLU_C': 6, 'GLU_O': 8, 'GLU_CB': 8, 'GLU_CG': 8, 'GLU_CD': 6, 'GLU_OE1': 8, 'GLU_OE2': 9, 'GLU_OXT': 9,
             'LYS_N': 8, 'LYS_CA': 7, 'LYS_C': 6, 'LYS_O': 8, 'LYS_CB': 8, 'LYS_CG': 8, 'LYS_CD': 8, 'LYS_CE': 8, 'LYS_NZ': 9, 'LYS_OXT': 9,
             'ARG_N': 8, 'ARG_CA': 7, 'ARG_C': 6, 'ARG_O': 8, 'ARG_CB': 8, 'ARG_CG': 8, 'ARG_CD': 8, 'ARG_NE': 8, 'ARG_CZ': 6, 'ARG_NH1': 8, 'ARG_NH2': 9, 'ARG_OXT': 9,
             'HIS_N': 8, 'HIS_CA': 7, 'HIS_C': 6, 'HIS_O': 8, 'HIS_CB': 8, 'HIS_CG': 6, 'HIS_ND1': 8, 'HIS_CD2': 7, 'HIS_CE1': 7, 'HIS_NE2': 8, 'HIS_OXT': 9}

atomType = {'GLY_N': 'N_single_bb', 'GLY_CA': 'C_single_bb', 'GLY_C': 'C_double_bb', 'GLY_O': 'O_double_bb', 'GLY_OXT': 'O_single',
            'ALA_N': 'N_single_bb', 'ALA_CA': 'C_single_bb', 'ALA_C': 'C_double_bb', 'ALA_O': 'O_double_bb', 'ALA_CB': 'C_single', 'ALA_OXT': 'O_single',
            'VAL_N': 'N_single_bb', 'VAL_CA': 'C_single_bb', 'VAL_C': 'C_double_bb', 'VAL_O': 'O_double_bb', 'VAL_CB': 'C_single', 'VAL_CG1': 'C_single', 'VAL_CG2': 'C_single', 'VAL_OXT': 'O_single',
            'LEU_N': 'N_single_bb', 'LEU_CA': 'C_single_bb', 'LEU_C': 'C_double_bb', 'LEU_O': 'O_double_bb', 'LEU_CB': 'C_single', 'LEU_CG': 'C_single', 'LEU_CD1': 'C_single', 'LEU_CD2': 'C_single', 'LEU_OXT': 'O_single',
            'ILE_N': 'N_single_bb', 'ILE_CA': 'C_single_bb', 'ILE_C': 'C_double_bb', 'ILE_O': 'O_double_bb', 'ILE_CB': 'C_single', 'ILE_CG1': 'C_single', 'ILE_CG2': 'C_single', 'ILE_CD1': 'C_single', 'ILE_OXT': 'O_single',
            'MET_N': 'N_single_bb', 'MET_CA': 'C_single_bb', 'MET_C': 'C_double_bb', 'MET_O': 'O_double_bb', 'MET_CB': 'C_single', 'MET_CG': 'C_single', 'MET_SD': 'S_single', 'MET_CE': 'C_single', 'MET_OXT': 'O_single',
            'PHE_N': 'N_single_bb', 'PHE_CA': 'C_single_bb', 'PHE_C': 'C_double_bb', 'PHE_O': 'O_double_bb', 'PHE_CB': 'C_single', 'PHE_CG': 'C_intermediate', 'PHE_CD1': 'C_intermediate', 'PHE_CD2': 'C_intermediate', 'PHE_CE1': 'C_intermediate', 'PHE_CE2': 'C_intermediate', 'PHE_CZ': 'C_intermediate', 'PHE_OXT': 'O_single',
            'TRP_N': 'N_single_bb', 'TRP_CA': 'C_single_bb', 'TRP_C': 'C_double_bb', 'TRP_O': 'O_double_bb', 'TRP_CB': 'C_single', 'TRP_CG': 'C_double', 'TRP_CD1': 'C_double', 'TRP_CD2': 'C_intermediate', 'TRP_NE1': 'N_single', 'TRP_CE2': 'C_intermediate', 'TRP_CE3': 'C_intermediate', 'TRP_CZ2': 'C_intermediate', 'TRP_CZ3': 'C_intermediate', 'TRP_CH2': 'C_intermediate', 'TRP_OXT': 'O_single',
            'PRO_N': 'N_single_bb', 'PRO_CA': 'C_single_bb', 'PRO_C': 'C_double_bb', 'PRO_O': 'O_double_bb', 'PRO_CB': 'C_single', 'PRO_CG': 'C_single', 'PRO_CD': 'C_single', 'PRO_OXT': 'O_single',
            'SER_N': 'N_single_bb', 'SER_CA': 'C_single_bb', 'SER_C': 'C_double_bb', 'SER_O': 'O_double_bb', 'SER_CB': 'C_single', 'SER_OG': 'O_single', 'SER_OXT': 'O_single',
            'THR_N': 'N_single_bb', 'THR_CA': 'C_single_bb', 'THR_C': 'C_double_bb', 'THR_O': 'O_double_bb', 'THR_CB': 'C_single', 'THR_OG1': 'O_single', 'THR_CG2': 'C_single', 'THR_OXT': 'O_single',
            'CYS_N': 'N_single_bb', 'CYS_CA': 'C_single_bb', 'CYS_C': 'C_double_bb', 'CYS_O': 'O_double_bb', 'CYS_CB': 'C_single', 'CYS_SG': 'S_single', 'CYS_OXT': 'O_single',
            'TYR_N': 'N_single_bb', 'TYR_CA': 'C_single_bb', 'TYR_C': 'C_double_bb', 'TYR_O': 'O_double_bb', 'TYR_CB': 'C_single', 'TYR_CG': 'C_intermediate', 'TYR_CD1': 'C_intermediate', 'TYR_CD2': 'C_intermediate', 'TYR_CE1': 'C_intermediate', 'TYR_CE2': 'C_intermediate', 'TYR_CZ': 'C_intermediate', 'TYR_OH': 'O_single', 'TYR_OXT': 'O_single',
            'ASN_N': 'N_single_bb', 'ASN_CA': 'C_single_bb', 'ASN_C': 'C_double_bb', 'ASN_O': 'O_double_bb', 'ASN_CB': 'C_single', 'ASN_CG': 'C_double', 'ASN_OD1': 'O_double', 'ASN_ND2': 'N_single', 'ASN_OXT': 'O_single',
            'GLN_N': 'N_single_bb', 'GLN_CA': 'C_single_bb', 'GLN_C': 'C_double_bb', 'GLN_O': 'O_double_bb', 'GLN_CB': 'C_single', 'GLN_CG': 'C_single', 'GLN_CD': 'C_double', 'GLN_OE1': 'O_double', 'GLN_NE2': 'N_single', 'GLN_OXT': 'O_single',
            'ASP_N': 'N_single_bb', 'ASP_CA': 'C_single_bb', 'ASP_C': 'C_double_bb', 'ASP_O': 'O_double_bb', 'ASP_CB': 'C_single', 'ASP_CG': 'C_double', 'ASP_OD1': 'O_intermediate', 'ASP_OD2': 'O_intermediate', 'ASP_OXT': 'O_single',
            'GLU_N': 'N_single_bb', 'GLU_CA': 'C_single_bb', 'GLU_C': 'C_double_bb', 'GLU_O': 'O_double_bb', 'GLU_CB': 'C_single', 'GLU_CG': 'C_single', 'GLU_CD': 'C_double', 'GLU_OE1': 'O_intermediate', 'GLU_OE2': 'O_intermediate', 'GLU_OXT': 'O_single',
            'LYS_N': 'N_single_bb', 'LYS_CA': 'C_single_bb', 'LYS_C': 'C_double_bb', 'LYS_O': 'O_double_bb', 'LYS_CB': 'C_single', 'LYS_CG': 'C_single', 'LYS_CD': 'C_single', 'LYS_CE': 'C_single', 'LYS_NZ': 'N_single', 'LYS_OXT': 'O_single',
            'ARG_N': 'N_single_bb', 'ARG_CA': 'C_single_bb', 'ARG_C': 'C_double_bb', 'ARG_O': 'O_double_bb', 'ARG_CB': 'C_single', 'ARG_CG': 'C_single', 'ARG_CD': 'C_single', 'ARG_NE': 'N_single', 'ARG_CZ': 'C_double', 'ARG_NH1': 'N_double', 'ARG_NH2': 'N_single', 'ARG_OXT': 'O_single',
            'HIS_N': 'N_single_bb', 'HIS_CA': 'C_single_bb', 'HIS_C': 'C_double_bb', 'HIS_O': 'O_double_bb', 'HIS_CB': 'C_single', 'HIS_CG': 'C_double', 'HIS_ND1': 'N_single', 'HIS_CD2': 'C_double', 'HIS_CE1': 'C_double', 'HIS_NE2': 'N_double', 'HIS_OXT': 'O_single'}

## Data from https://arxiv.org/pdf/0804.2488.pdf
radii = {'C_single': 0.77, 'C_double': 0.67, 'C_intermediate': 0.72, 'C_single_bb': 0.77, 'C_double_bb': 0.67,
         'O_single': 0.67, 'O_double': 0.60, 'O_intermediate': 0.64, 'O_double_bb': 0.60,
         'N_single': 0.70, 'N_double': 0.62, 'N_intermediate': 0.66, 'N_single_bb': 0.70,
         'S_single': 1.04}

totalElectrons = {}
for atom, num in electrons.items():
    if atom[-3:] == 'OXT': continue

    if atom[0:3] in totalElectrons.keys():
        totalElectrons[atom[0:3]] += num
    else:
        totalElectrons[atom[0:3]] = num

pdbidfile = sys.argv[1]
pdbids = []
with open(pdbidfile, "r") as fileHandle:
    for pdbid in fileHandle:
        pdbids.append(pdbid.split(" ; ")[0])
fileHandle.close()

fileHandle = open(sys.argv[2], 'w')
#print(*["pdbid", "resolution", "spaceGroup", "chainMean", "chainMedian", "chainLogMean", "chainLogMedian", "resMean",
#        "resMedian", "resLogMean", "resLogMedian", "atomMean", "atomMedian", "atomLogMean", "atomLogMedian"], sep=", ",
#      file=fileHandle)

for pdbid in pdbids:
    try:
        pdbid = pdbid.lower()
        densityObj = ccp4.readFromPDBID(pdbid)
        densityCutoff = sigma = np.mean(densityObj.densityArray) + 1.5 * np.std(densityObj.densityArray)

        pdbfile = './pdb/pdb' + pdbid + '.ent'
        if os.path.isfile(pdbfile):
            pass
        else:
            pdbl = pdb.PDBList()
            pdbl.retrieve_pdb_file(pdbid, pdir='./pdb/')

        # Bio Python parser
        parser = pdb.PDBParser()
        structure = parser.get_structure(pdbid, pdbfile)

        ## my own parser
        pdbObj = myPDB.readPDBfile(pdbfile)
        program = pdbObj.header.program
        spaceGroup = pdbObj.header.spaceGroup

        '''
        ## crystal contacts by Michael
        mmcifFile = '/mlab/data/databases/PDB/mmcif/2015_09_04/' + pdbid[1:3] + '/' + pdbid + '.cif.gz'
        if os.path.isfile(mmcifFile):
            pass
        else:
            print("No", pdbid, "mmcif file exist" )
            continue

        cc = crystalContacts.get_contact_atoms(mmcifFile)
        contactAtoms = [str(atom.parent.id[1]) + '_' + atom.name for atom in cc]

        atomNum = sum(1 for x in structure.get_atoms())
        contactNum = sum(1 for x in contactAtoms)
        '''
    except:
        continue

    # if "organism_scientific" not in structure.header["source"]["1"].keys():
    #    continue

    # atoms = structure.get_atoms()
    # atoms = list(structure.get_atoms())
    # atoms[0].get_coord()

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

            blobs = densityObj.findAberrantBlobs(atom.coord, radii[atomType[resAtom]], densityCutoff)
            if len(blobs) == 0:
                continue

            #atomList.append(str(residue.id[1]) + ', ' + resAtom + ', ' + atomType[resAtom] + ', ' + str(blobs[0].totalDensity / electrons[resAtom]) + ', ' + str(atom.occupancy) + ', ' + str(len(blobs[0].crsList)))
            atomList.append([residue.id[1], resAtom, atomType[resAtom], blobs[0].totalDensity / electrons[resAtom], len(blobs[0].crsList)])
            atomAvgDensity.append(blobs[0].totalDensity / electrons[resAtom])


            for blob in blobs:
                for cloud in resClouds:
                    if cloud.testOverlap(blob):
                        atoms = cloud.atoms
                        cloud.merge(blob, densityObj.density)
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
            resAvgDensity.append(
                cloud.totalDensity / sum([electrons[atom.parent.resname + '_' + atom.name] for atom in cloud.atoms]))

        for blob in resClouds:
            for cloud in chainClouds:
                if cloud.testOverlap(blob):
                    atoms = cloud.atoms
                    cloud.merge(blob, densityObj.density)
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

    """
    # print(cloud.centroid, cloud.volume, cloud.totalDensity, totalElectron, [x.serial_number for x in cloud.atoms], file=fileHandle)

    # for resname, blobList in resDict.items():
    #    print(resname, len(blobList), [x.totalDensity for x in blobList], sum([x.totalDensity for x in blobList])/len(blobList)/totalElectrons[resname], file=fileHandle)
    # print(*resList, sep="; ")

    chainLogMean = 10 ** np.mean([np.log10(x) for x in chainAvgDensity])
    chainLogMedian = 10 ** np.median([np.log10(x) for x in chainAvgDensity])
    chainMean = np.mean(chainAvgDensity)
    chainMedian = np.median(chainAvgDensity)

    resLogMean = 10 ** np.mean([np.log10(x) for x in resAvgDensity])
    resLogMedian = 10 ** np.median([np.log10(x) for x in resAvgDensity])
    resMean = np.mean(resAvgDensity)
    resMedian = np.median(resAvgDensity)

    atomLogMean = 10 ** np.mean([np.log10(x) for x in atomAvgDensity])
    atomLogMedian = 10 ** np.median([np.log10(x) for x in atomAvgDensity])
    atomMean = np.mean(atomAvgDensity)
    atomMedian = np.median(atomAvgDensity)

    print(*[pdbid, structure.header["resolution"], spaceGroup, chainMean, chainMedian, chainLogMean, chainLogMedian,
            resMean, resMedian, resLogMean, resLogMedian, atomMean, atomMedian, atomLogMean, atomLogMedian], sep=", ",
          file=fileHandle)
    """

    # reduce atom type to element and single/double/intermediate
    def atomtype(x):
        aa = x.split('_')
        return aa[0] + '_' + aa[1]

    # normalize the density by median volumn of given atom type
    def normVolumn(row):
        return float(row['density']) / float(row['volumn']) * float(medians['volumn'][row['atomType']])

    atoms = pandas.DataFrame(atomList, columns=['resID', 'resAtom', 'atomType', 'density', 'volumn'])
    atoms['atomType'] = atoms['atomType'].apply(atomtype)
    medians = atoms.groupby(['atomType']).median()
    atoms['adjDensity'] = atoms.apply(lambda row: normVolumn(row), axis=1)
    medians = atoms.groupby(['atomType']).median()
    medianAdjDen = medians.sort_index()['adjDensity']

    print(*medianAdjDen, sep=", ", file=fileHandle)

fileHandle.close()
