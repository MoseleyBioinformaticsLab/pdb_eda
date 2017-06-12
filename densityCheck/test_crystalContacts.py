#!/usr/bin/python3

import Bio.PDB as pdb
import crystalContacts

pdbid = "1g3v"
pdbfile = './pdb/pdb' + pdbid.lower() + '.ent'

parser = pdb.PDBParser()
structure = parser.get_structure(pdbid, pdbfile)

atoms = crystalContacts.get_contact_atoms(structure)
print(atoms[0])