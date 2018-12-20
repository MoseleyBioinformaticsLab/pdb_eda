import os
import sys
import math
import time

import numpy as np
from scipy import stats
import multiprocessing
import datetime
import tempfile

from . import densityAnalysis

## Final set of radii derived from optimization
radiiDefault = {'C_single': 0.84, 'C_double': 0.67, 'C_intermediate': 0.72, 'C_single_bb': 0.72, 'C_double_bb': 0.61, 
                'O_single': 0.80, 'O_double': 0.77, 'O_intermediate': 0.82, 'O_double_bb': 0.71, 
                'N_single': 0.95, 'N_intermediate': 0.77, 'N_single_bb': 0.7, 
                'S_single': 0.75}

slopesDefault = {'C_single': -0.33373359301010996, 'C_double': -0.79742538209281033, 'C_intermediate': -0.46605044936397311, 'C_single_bb': -0.44626029983492604, 'C_double_bb': -0.53313491315106321, 
                'O_single': -0.61612691527685515, 'O_double': -0.69081386892018048, 'O_intermediate': -0.64900152439990022, 'O_double_bb': -0.59780395867126568, 
                'N_single': -0.50069328952200343, 'N_intermediate': -0.5791543104296577, 'N_single_bb': -0.4905830761073553, 
                'S_single': -0.82192528851026947}

def processFunction(pdbid, radii, slopes):
    analyser = densityAnalysis.fromPDBid(pdbid)

    if not analyser:
            return 0 

    analyser.aggregateCloud(radii, slopes)
    if not analyser.chainMedian:
        return 0

    diffs = []
    stats = []
    for atomType in sorted(radiiDefault):
        diff = (analyser.medians.loc[atomType]['correctedDensity'] - analyser.chainMedian) / analyser.chainMedian if atomType in analyser.medians.index else 0
        diffs.append(diff)
    stats = [analyser.pdbid, analyser.chainMedian, analyser.pdbObj.header.resolution, analyser.pdbObj.header.spaceGroup]
    #return diffs, stats

    globalTempFile.write("%s\n" % ', '.join([str(i) for i in stats + diffs]))
    return globalTempFile.name

def openTemporaryFile(temp):
    dirname = os.getcwd()
    global globalTempFile
    globalTempFile = tempfile.NamedTemporaryFile(mode='w', dir=dirname, prefix="tempPDB_", delete=False)
    time.sleep(0.01) # sleep 0.01 seconds to prevent the same worker from calling this twice.
    return globalTempFile.name

def closeTemporaryFile(temp):
    filename = globalTempFile.name
    globalTempFile.close()
    # Sleep 0.01 seconds to prevent the same worker from calling this twice.
    # May need to increase this to 1 second or even 10 seconds to make this work.
    time.sleep(0.01)
    return filename


def main(pdbidfile, resultname, mode):
    '''
    Mode: 
        0 - use optimized radii
        1 - use original radii
    '''
    pdbids = []
    with open(pdbidfile, "r") as fileHandleIn:
        for pdbid in fileHandleIn:
            pdbids.append(pdbid[0:4])

    radii = radiiDefault
    slopes = slopesDefault
    if int(mode) == 0:
        radii = {}
        slopes = {}

    with multiprocessing.Pool() as pool:
        result_filenames = pool.map(openTemporaryFile, (1 for i in range(multiprocessing.cpu_count())))
        results = pool.starmap(processFunction, ((pdbid, radii, slopes) for pdbid in pdbids))
        closed_filenames = pool.map(closeTemporaryFile, (1 for i in range(multiprocessing.cpu_count())))

    #fileHandle = open(resultname, 'w')
    #for result in results:
    #    if result:
    #        print(*result[1], *result[0], sep=", ", file=fileHandle)

    unclosed_filenames = set(result_filenames) - set(closed_filenames)
    if unclosed_filenames: # check whether all result files were closed.  Increase sleep time in closeTemporary File to prevent this.
        print("Unclosed Files: ", unclosed_filenames)

    with open(resultname, "w") as outfile:
        for f in result_filenames:
            if f:
                with open(f, "r") as infile:
                    outfile.write(infile.read())
                os.remove(f) # remove the file

    '''
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

def singleCalc(pdbid, resultDir, mode):
    '''
    Mode: 0 - use optimized radii
          1 - use original radii
    '''
    radii = radiiDefault
    slopes = slopesDefault
    if int(mode) == 0:
        radii = {}
        slopes = {}

    analyser = densityAnalysis.fromPDBid(pdbid)
    analyser.aggregateCloud(radii, slopes, atomL=True, residueL=True, chainL=True)

    analyser.atomList.to_csv(resultDir + pdbid + '.atom.txt', sep=',')
    
    fileHandle = open(resultDir + pdbid + '.residue.txt', 'w')
    for item in analyser.residueList:
        print(', '.join(map(str, item)), file=fileHandle)
    fileHandle.close()
    
    fileHandle = open(resultDir + pdbid + ".chain.txt", 'w')
    #print(*analyser.chainList, sep='\n', file=fileHandle)
    for item in analyser.chainList:
        print(', '.join(map(str, item)), file=fileHandle)
    fileHandle.close()

if __name__ == '__main__':
    """
    runMode:
        s - single structure
        m - multiple structures
    radiiMode: 
        0 - use optimized radii
        1 - use original radii
    """

    if len(sys.argv) == 5:
        _, runMode, radiiMode, filename, resultname = sys.argv
        if runMode == 's' or runMode == 'single':
            singleCalc(filename, resultname, radiiMode)
        elif runMode == 'm' or runMode == 'multiple':
            main(filename, resultname, radiiMode)
        else:
            print("Wrong options")
    else:
        print("Wrong number of arguments!")



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
