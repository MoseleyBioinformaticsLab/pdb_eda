import pdb

pdbFile = "./pdb/pdb1g3v.ent"
pdbObj = pdb.readPDBfile(pdbFile)

print(pdbObj.header.pdbid, pdbObj.header.program, pdbObj.header.spaceGroup, pdbObj.header.date, pdbObj.header.method, pdbObj.header.resolution, pdbObj.header.rValue, pdbObj.header.rFree, "end")