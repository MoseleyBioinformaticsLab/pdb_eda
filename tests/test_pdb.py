from ..pdb_eda import pdb

pdbFile = "./pdb/pdb1g3v.ent"
pdbObj = pdb.readPDBfile(pdbFile)

print(pdbObj.header.pdbid)