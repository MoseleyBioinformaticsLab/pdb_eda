from pdb_eda import pdbParser
import urllib.request

pdburlPrefix = "http://ftp.rcsb.org/pub/pdb/data/structures/all/pdb/"
pdbfile = 'pdb1g3v.ent.gz'
url = pdburlPrefix + pdbfile
urllib.request.urlretrieve(url, pdbfile)

with gzip.open(pdbfile, 'rt') as gzipFile:
    pdbObj = pdbParser.readPDBfile(gzipFile)

print(pdbObj.header.pdbid)