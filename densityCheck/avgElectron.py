# !/usr/bin/python3

import Bio.PDB as pdb
import ccp4

pdbid = '1cbs'
densityObj = ccp4.readFromPDBID(pdbid)

parser = pdb.PDBParser()
structure = parser.get_structure(pdbid, pdbid+'.pdb')

atoms = structure.get_atoms()
#atoms = list(structure.get_atoms())
#atoms[0].get_coord()

radius = 0.8

clouds = []
for atom in structure.get_atoms():
    blobs = densityObj.findAberrantBlobs(atom.coord, radius)

    for blob in blobs:
        for cloud in clouds:
            if cloud.testOverlap(blob):
                atoms = cloud.atoms
                cloud.merge(blob, densityObj.density)
                cloud.atoms = atoms + [atom]
                break

        else:
            blob.atoms = [atom]
            clouds.append(blob)


elementElectron = {'C': 6, 'N': 7, 'O': 8, 'P': 15, 'S': 16}

for cloud in clouds:
    cloud.totalElectron = sum([elementElectron[atom.element] for atom in cloud.atoms])

fileHandle = open('output.txt', 'w')
for cloud in clouds:
    print(cloud.centroid, cloud.volume, cloud.totalDensity, cloud.totalElectron, [x.serial_number for x in cloud.atoms], file=fileHandle)
fileHandle.close()

