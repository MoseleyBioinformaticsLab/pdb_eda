from TEMPy.StructureParser import PDBParser
structure_instance=PDBParser.fetch_PDB('1cbs','1cbs.pdb',hetatm=True,water=False)

from TEMPy.MapParser import MapParser
target_map=MapParser.readMRC('1cbs.ccp4')
target_map.x_origin()
target_map.getMap()
aa = target_map.getMap()
aa.shape
aa = aa.tolist()
hist = plt.hist(aa, 200)
plt.show()