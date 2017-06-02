from ..densityCheck import pdb

pdbFile = "./pdb/pdb1g3v.ent"
pdbObj = pdb.readPDBfile(pdbFile)

print(pdbObj.header.pdbid)