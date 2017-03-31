# !/usr/bin/python3

import ccp4
import matplotlib.pyplot as plt
import numpy as np

#pdbid = '1cbs' # 2Fo - Fc
#density1 = ccp4.readFromPDBID(pdbid)

pdbid = '1cbs_diff' # Fo - Fc
densityObj = ccp4.readFromPDBID(pdbid)

testPos = [6.1, 21.3, 20.7]  # green, big positive
testPos1 = [34.3, 25.3, 22.8]  # red, big negative
print('density: ', densityObj.getPointDensityFromXyz(testPos1))

cutoff = -3 * densityObj.header.rmsd
blub = densityObj.findAberrantBlubs(testPos1, 5, cutoff)
print('aberrant blubs: ', blub)
print('test position: ', testPos1)
print('unit volume: ', densityObj.header.unitVolume)
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


