#!/usr/bin/python3

import Bio.PDB as pdb
import crystalContacts

pdbid = "1g3v"
#pdbfile = './pdb/pdb' + pdbid.lower() + '.ent'
pdbfile = '/mlab/data/databases/PDB/mmcif/2015_09_04/g3/1g3z.cif.gz'
mmcifFile = '/mlab/data/databases/PDB/mmcif/2015_09_04/' + pdbid.lower()[1:3] + '/' + pdbid + '.cif.gz'

atoms = crystalContacts.get_contact_atoms(mmcifFile)
atom = atoms[100]

print(len(atoms))
print(atom.parent.resname + '_' + atom.name, atom.parent.id[1])
print(atom.parent.__dict__)
print(atom.__dict__)
