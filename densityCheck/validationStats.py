# !/usr/bin/python3
"""
validationStats.py
    Calculate real space correlation coefficients (RSCC) and real space R value (RSR) via different methods.
"""
import os.path
#import requests
import numpy as np
from scipy.stats.stats import pearsonr


class validationStats(object):
    def __init__(self, pdbid):
        self.pdbid = pdbid


    def getStats(self, structure, fo, fc, sigma):
        """
        RETURNS the rscc values from the Fo and Fc density maps
        """
        resolution = structure.header['resolution']
        radius = 0.7
        if 0.6 <= resolution <= 3:
            radius = (resolution - 0.6) / 3 + 0.7
        elif resolution > 3:
            radius = resolution * 0.5

        statsList = []
        for residue in structure.get_residues():
            if residue.id[0] != ' ':
                continue

            crsLists = []
            bfactor = occupancy = 0
            for atom in residue.child_list:
                crsList = fo.getSphereCrsFromXyz(atom.coord, radius, sigma)
                crsLists = crsLists + [i for i in crsList if i not in crsLists]

                bfactor = bfactor + atom.get_bfactor() * atom.get_occupancy()
                occupancy = occupancy + atom.get_occupancy()

            foDensity = [fo.getPointDensityFromCrs(i) for i in crsLists]
            fcDensity = [fc.getPointDensityFromCrs(i) for i in crsLists]

            rscc = pearsonr(foDensity, fcDensity)[0]
            rsr = sum(abs(np.array(foDensity) - np.array(fcDensity))) / sum(abs(np.array(foDensity) + np.array(fcDensity)))

            # residue num, residue name, rscc, rsr, occupancy-weighted average B factor, number of involving grid points
            statsList.append(", ".join([residue.parent.id, str(residue.id[1]), residue.resname, str(rscc), str(rsr), str(bfactor / occupancy), str(len(crsLists))]))

        return statsList


    def getEDSstats(self):
        """
        RETURNS the rscc values from Uppsala Electron Density Server (EDS)
        """
        statsFilePath = '/mlab/project/metal/rscc/stats/' + self.pdbid + '_stat.lis'
        statsFileUrl = str('http://eds.bmc.uu.se/eds/dfs/' + self.pdbid[1:-1] + '/' + self.pdbid + '/' + self.pdbid + 'stat.lis').lower
        response = requests.get(statsFileUrl)
        if os.path.isfile(statsFilePath):
            statsFile = open(statsFilePath)
            lines = statsFile.read().split("\n")
        elif response.status_code == 200:
            lines = response.text.split('\n')
        else:
            return 0

        statsList = []
        for line in lines:
            if line[0] == '!' or line[21:26] != ' ':
                continue
            else:
                rscc = line[21:26]
                rsr = line[27:32]
                # chain id, residue num, residue name, rscc, occupancy-weighted average B factor, number of grid points
                statsList.append(", ".join([line[12:13], line[15:17], line[8:11], rscc, rsr, line[35:40], line[80:83]]))

        return statsList
