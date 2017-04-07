# !/usr/bin/python3

import ccp4
import matplotlib.pyplot as plt
import numpy as np

# pdbid = '1cbs' # 2Fo - Fc
# density1 = ccp4.readFromPDBID(pdbid)

pdbid = '1cbs_diff'  # Fo - Fc
densityObj = ccp4.readFromPDBID(pdbid)
print("origin: ", densityObj.origin)

testPos = [6.1, 21.3, 20.7]  # green, big positive
testPos1 = [34.3, 25.3, 22.8]  # red, big negative
print('density: ', densityObj.getPointDensityFromXyz(testPos1))

cutoff = 3 * densityObj.header.rmsd

densityObj.density[:20, :20, :20] = 0
centerXYZ = densityObj.header.crs2xyzCoord([10, 10, 10])

# Set up green blob
densityObj.density[11:13, 11:13, 11:13] = 1  # (x,y,z), Volume 2*2*2=8
densityObj.density[7:10, 11:14, 7:10] = 1  # (-x,y,-z), Volume 3*3*3=27

# Set up red blob
densityObj.density[11:13, 7:9, 7:9] = -1  # (-x,-y,z), Volume 2*2*2=8
densityObj.density[7:10, 7:10, 11:14] = -1  # (x,-y,-z), Volume 3*3*3=27

blub = densityObj.findAberrantBlobs(centerXYZ, 5, -cutoff)

print('unit volume: ', densityObj.header.unitVolume)
print('aberrant blubs: ', blub)

print(densityObj.header.crs2xyzCoord([12, 8, 8]))

'''
fo = np.array(density1) - np.array(density2)
fc = np.array(density1) - np.array(density2) * 2

bins = np.linspace(-1.5, 3, 500)

plt.hist(fo, bins, alpha=0.5, label='fo')
plt.hist(fc, bins, alpha=0.5, label='fc')
plt.legend(loc='upper right')
plt.show()
'''

'''
hist = plt.hist(densities, 200)
idx = hist[0].tolist().index(max(hist[0]))
print(hist[1])
print(idx, hist[0][idx], hist[1][idx], hist[1][idx+1], hist[0][0], hist[1][0])
plt.show()
'''

'''
#fileName = "../1cbs_diff.ccp4"
fileName = "../1cbs.ccp4"
ccp4.read(fileName)
'''


