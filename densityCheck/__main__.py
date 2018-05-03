import sys
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

from . import densityAnalysis


def main(pdbidfile):
    pdbids = []
    with open(pdbidfile, "r") as fileHandleIn:
        for pdbid in fileHandleIn:
            pdbids.append(pdbid.split(" ; ")[0])

    fileHandle = open("results/densityCutoffs.2.txt", 'w')
    for pdbid in pdbids:
        analyser = densityAnalysis.fromPDBid(pdbid)

        if not analyser:
            continue

        print(*[pdbid, analyser.densityObj.densityCutoffFromHeader, analyser.densityObj.densityCutoff, analyser.densityObj.densityCutoffFromLeftSide, analyser.densityObj.densityCutoffFromLeftSide2, analyser.densityObj.densityCutoffFromLeftSide25, analyser.densityObj.densityCutoffFromLeftSide3], sep=', ', file=fileHandle)
    fileHandle.close()

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


        # diffDenStats = analyser.calcAtomBlobDists()
        # if math.isnan(analyser.chainMedian) or analyser.chainMedian is None:
        #     print('flag')
        #     continue
        # dists = [i[0] for i in diffDenStats]
        #
        # '''
        # # in normal space
        # kernel = stats.gaussian_kde(dists)
        # x = np.linspace(min(dists), max(dists), 200)
        # mode = x[np.argmax(kernel(x))]
        # leftside = [i for i in dists if i < mode]
        # dev = np.sqrt(sum([(i - mode) ** 2 for i in leftside]) / len(leftside))
        # cutoff = mode + dev * 2
        # cutoff1 = mode + dev * 3
        # '''
        #
        # # in log space
        # logdists = np.log(dists)
        # kernel = stats.gaussian_kde(logdists)
        # x = np.linspace(min(logdists), max(logdists), 200)
        # logmode = x[np.argmax(kernel(x))]
        # logleftside = [i for i in logdists if i < logmode]
        # logdev = np.sqrt(sum([(i - logmode) ** 2 for i in logleftside]) / len(logleftside))
        # logcutoff = logmode + logdev * 2
        # #logcutoff1 = logmode + logdev * 3
        #
        # isolatedcutoff = 5
        # if np.exp(logcutoff) < 5:
        #     isolatedcutoff = np.exp(logcutoff)
        # electrons = [i[2] for i in diffDenStats]
        # nonIsolatedElectrons = [e for i, e in enumerate(electrons) if dists[i] <= isolatedcutoff]
        # isolatedElectrons = [e for i, e in enumerate(electrons) if dists[i] > isolatedcutoff and diffDenStats[i][1] == '-']
        #
        # if len(nonIsolatedElectrons)  == 0 or len(isolatedElectrons) == 0:
        #     continue
        #
        # f, axarr = plt.subplots(2, sharex=True)
        # axarr[0].hist(nonIsolatedElectrons, 100)
        # axarr[0].set_title('non-isolated electrons')
        # axarr[1].hist(isolatedElectrons, 50)
        # axarr[0].set_title('isolated red blob electrons')
        # plt.savefig('../electrons/' + pdbid + '.png')
        # plt.close()


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

if __name__ == '__main__':
    _, filename = sys.argv

    main(filename)

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
