import sys
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

from . import densityAnalysis


def main(pdbidfile, resultname, atomType, radius):
    pdbids = []
    with open(pdbidfile, "r") as fileHandleIn:
        for pdbid in fileHandleIn:
            pdbids.append(pdbid[0:4])

    fileHandle = open(resultname, 'w')
    n = 0
    diff = []
    for pdbid in pdbids:
        analyser = densityAnalysis.fromPDBid(pdbid)

        if not analyser:
            continue

        analyser.aggregateCloud(atomType, radius)
        if not analyser.chainMedian:
            continue
        if n == 0: 
            n = 1
            print("pdbid", "chainMedian", *list(analyser.medians.index), *list(analyser.medians.index), sep=', ', file=fileHandle)
        print(pdbid, analyser.chainMedian, *analyser.medians['correctedDensity'], *analyser.medians['slopes'], sep=", ", file=fileHandle) 
        if atomType in analyser.medians.index:
            diff.append((analyser.medians.loc[atomType]['correctedDensity'] - analyser.chainMedian) / analyser.chainMedian)  ## for radii optimization

        '''
        for item in analyser.atomList:
            print(', '.join(map(str, item + [analyser.chainMedian])), file=fileHandle)
        #print(pdbid + ', ' + str(analyser.chainMedian), file=fileHandle)
        fileHandle.close()

        if math.isnan(analyser.chainMedian) or analyser.chainMedian is None:
            print(pdbid + ' has no aggregated chain clouds.')
            continue
        dists = [i[0] for i in diffDenStats]
        electrons = [i[2] for i in diffDenStats]
        npoints = [i[3] for i in diffDenStats]

        # in log space
        logdists = np.log(dists)
        kernel = stats.gaussian_kde(logdists)
        x = np.linspace(min(logdists), max(logdists), 200, dtype=np.float64)
        logmode = x[np.argmax(kernel(x))]
        logleftside = [i for i in logdists if i < logmode]
        logdev = np.sqrt(sum([(i - logmode) ** 2 for i in logleftside]) / len(logleftside))
        logcutoff = logmode + logdev * 2

        isolatedEcutoff = 5
        if np.exp(logcutoff) < 5:
            isolatedEcutoff = np.exp(logcutoff)
        nonIsolatedElectrons = [e for i, e in enumerate(electrons) if dists[i] <= isolatedEcutoff]
        isolatedElectrons = [e for i, e in enumerate(electrons) if dists[i] > isolatedEcutoff and diffDenStats[i][1] == '-']
        randomElectrons = [e for i, e in enumerate(electrons) if dists[i] > isolatedEcutoff and diffDenStats[i][1] == '-' and npoints[i] == 1]
        intrinsicElectrons = [e for i, e in enumerate(electrons) if dists[i] > isolatedEcutoff and diffDenStats[i][1] == '-' and npoints[i] > 1]


        # if len(nonIsolatedElectrons)  == 0 or len(isolatedElectrons) == 0:
        #     continue

        plt.figure(figsize=(15, 10), dpi=80)
        ax1 = plt.subplot(421)
        ax2 = plt.subplot(423, sharex=ax1)
        ax3 = plt.subplot(425, sharex=ax1)
        ax4 = plt.subplot(427, sharex=ax1)
        ax5 = plt.subplot(122)

        ax1.hist(nonIsolatedElectrons, 100)
        ax1.set_title('non-isolated electrons')
        ax2.hist(isolatedElectrons, 50)
        ax2.set_title('isolated red blob electrons')
        ax3.hist(randomElectrons, 50)
        ax3.set_title('random error electrons')
        ax4.hist(intrinsicElectrons, 50)
        ax4.set_title('intrinsicElectrons error electrons')
        ax5.hist(dists, 200)
        ax5.axvline(x=np.exp(logmode), color='black')
        ax5.axvline(x=np.exp(logcutoff), color='green')
        ax5.axvline(x=5, color='yellow')
        ax5.set_title('distance to the closest atom')

        plt.savefig('../error_properties/' + pdbid + '.png')
        plt.close()


        # plot cutoff
        # plt.hist(dists, 200)
        # plt.axvline(x=mode, color='black')
        # plt.axvline(x=cutoff, color='red')
        # plt.axvline(x=cutoff1, color='red')
        # plt.axvline(x=np.exp(logcutoff), color='green')
        # plt.axvline(x=np.exp(logcutoff1), color='green')
        # plt.axvline(x=5, color='yellow')
        # plt.savefig('../isolated_blob_cutoff/' + pdbid + '.png')
        # plt.close()

        # # testing different density cutoff
        # print(*[pdbid, analyser.densityObj.densityCutoffFromHeader, analyser.densityObj.densityCutoff, analyser.densityObj.densityCutoffFromLeftSide, analyser.densityObj.densityCutoffFromLeftSide2, analyser.densityObj.densityCutoffFromLeftSide25, analyser.densityObj.densityCutoffFromLeftSide3], sep=', ', file=fileHandle)

        # # validate the new aggregation method
        # analyser.aggregateCloud(atomL=True, residueL=True, chainL=True)
        #
        # fileHandle = open("results/aggregation/"+ pdbid + ".atomList.new.txt", 'w')
        # print(*analyser.atomList, sep='\n', file=fileHandle)
        # fileHandle.close()
        #
        # fileHandle = open("results/aggregation/" + pdbid + ".residueList.new.txt", 'w')
        # print(*analyser.residueList, sep='\n', file=fileHandle)
        # fileHandle.close()
        #
        # fileHandle = open("results/aggregation/" + pdbid + ".chainList.new.txt", 'w')
        # print(*analyser.chainList, sep='\n', file=fileHandle)
        # fileHandle.close()
    #fileHaddndle.close()
        '''

    print(np.nanmean(diff), np.nanmedian(diff), file=fileHandle) ## for radii optimization

if __name__ == '__main__':
    _, filename, resultname, atomType, radius = sys.argv

    main(filename, resultname, atomType, radius)

    #print(diffDenStats[111])
    #dists = [i[0] for i in diffDenStats]
    #plt.hist(dists, 200)




#suffix = sys.argv[2]
#fileHandle = open(sys.argv[2], 'w')
#radii[sys.argv[3]] = float(sys.argv[4])  # for radii optimization
#fileHandleB = open(sys.argv[3], 'w') #for b factor print out


# fileHandle = open("results/rscc." + pdbid + "." + suffix + ".txt", 'w')
# print(*rsccList, sep="\n", file=fileHandle)
# fileHandle.close()


# if n == 1:
# n = 0
# print("pdbid", "chainMedian", *sorted(radii.keys()), sep=', ', file=fileHandle)
# print("pdbid", "chainMedian", *sorted(radii.keys()), sep=', ', file=fileHandleB) ## for print out b factors
# print("pdbid", *sorted(atomTypeCount.keys()), sep=', ', file=fileHandle) ## for radii optimization, old

# print(pdbid, chainMedian, *medianAdjDen, sep=", ", file=fileHandle) ## for checking the medians of chain and all atom types
# print(pdbid, bfactorMedian, *bfactors, sep=", ", file=fileHandleB) ## for print out b factors
# print(pdbid, *[atomTypeCount[key] for key in sorted(atomTypeCount.keys())], sep=', ', file=fileHandle) ## for atom type composition

# diff.append((chainMedian - medianAdjDen[int(sys.argv[5])]) / chainMedian)  ## for radii optimization


#dists = [i[0] for i in diffMapStats]
#histogram = plt.hist(dists, 200)
#plt.savefig('../' + pdbid + '.png')
#plt.close()

#diffMapStats.sort(key=lambda x: x[2], reverse=True)  # sort by number of electron
#diffMapStats.sort(key=lambda x: x[4])  # sort by distance

#dists = [row[4] for row in diffMapStats]
#plt.hist(dists, bins=200)

#dists = [row[6] for row in diffMapStats]
#plt.hist(dists, bins=200)
#plt.savefig('../' + pdbid + '.png')
#plt.close()

#print(len([i for i in dists if i <= 2]))
#print(len(dists))

'''
dists = [i[1] for i in diffMapStats]
histogram = plt.hist(dists, bins=np.arange(min(dists), max(dists) + 0.02, 0.02))
plt.close()
mode = histogram[1][np.argmax(histogram[0][0:250]) + 1]  # maximum bin < 5A (5/0.02 = 250) distance
logmode = np.log(mode)

logdists = [np.log(i) for i in dists]
leftside = [i for i in logdists if i < logmode]
dev1 = np.sqrt(sum([(i - logmode) ** 2 for i in leftside]) / len(leftside))

bothside = [i for i in logdists if i < logmode + 2 * dev1]
cutoff = np.mean(bothside) + 2 * np.std(bothside)

plt.hist(logdists, 200)
plt.axvline(x=logmode, color='red')
plt.axvline(x=logmode + dev1, color='orange')
plt.axvline(x=logmode + 2 * dev1, color='yellow')
plt.axvline(x=cutoff, color='green')
plt.savefig('../atom-blue-distance/' + pdbid + '.log.png')
plt.close()

model = [i for i in dists if i < 10]
plt.hist(model, 200)
plt.axvline(x=mode, color='red')
plt.axvline(x=np.exp(logmode + dev1), color='orange')
plt.axvline(x=np.exp(logmode + 2 * dev1), color='yellow')
plt.axvline(x=np.exp(cutoff), color='green')
plt.savefig('../atom-blue-distance/' + pdbid + '.original.png')
plt.close()

model = [i for i in diffMapStats if i[1] < np.exp(cutoff)]

plt.hist([i[2] for i in model], bins=100)
plt.axvline(x=np.exp(cutoff), color='green')
plt.savefig('../min-red-blue-distance/' + pdbid + '.2.png')
plt.close()

plt.hist([i[3] for i in model], bins=100)
plt.axvline(x=np.exp(cutoff), color='green')
plt.savefig('../red-atom-avg-distance/' + pdbid + '.2.png')
plt.close()

plt.hist([i[4] for i in model], bins=100)
plt.axvline(x=np.exp(cutoff), color='green')
plt.savefig('../min-red-atom-distance/' + pdbid + '.2.png')
plt.close()
'''

#print(np.nanmean(diff), np.nanmedian(diff), file=fileHandle) ## for radii optimization

#fileHandle = open("results/cen." + pdbid + ".txt", 'w') ## form print out centroid-coordinates difference
#print(*atomList, sep='\n', file=fileHandle) ## for single atoms

#fileHandle.close()
#fileHandleB.close()

#print('Done! ', str(datetime.datetime.now()))
