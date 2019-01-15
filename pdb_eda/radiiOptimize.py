import sys
import math
import numpy as np
import multiprocessing
import datetime
from scipy import stats

from . import densityAnalysis

## Radii and slopes from an intial analysis on 100 structures
radiiDefault = {'C_single': 0.84, 'C_double': 0.66, 'C_intermediate': 0.70, 'C_single_bb': 0.69, 'C_double_bb': 0.61, 
                'O_single': 0.80, 'O_double': 0.77, 'O_intermediate': 0.82, 'O_double_bb': 0.71,
                'N_single': 0.90, 'N_intermediate': 0.75, 'N_single_bb': 0.69,
                'S_single': 0.75}

slopesDefault = {'C_double': -0.6538044, 'C_double_bb': -0.4626215, 'C_intermediate': -0.4494971, 'C_single': -0.3387809, 'C_single_bb': -0.3808402,
                'N_intermediate': -0.5541342, 'N_single': -0.4889789, 'N_single_bb': -0.5110914,
                'O_double': -0.7432083, 'O_double_bb': -0.6818212, 'O_intermediate': -0.7026091, 'O_single': -0.7070469,
                'S_single': -0.8644369}

def processFunction(pdbid, radii, slopes):
    analyser = densityAnalysis.fromPDBid(pdbid)

    if not analyser:
        return 0 

    analyser.aggregateCloud(radii, slopes)
    if not analyser.chainMedian:
        return 0

    diffs = []
    slopes = []
    for atomType in sorted(radiiDefault):
        #diff = (analyser.medians.loc[atomType]['correctedDensity'] - analyser.chainMedian) / analyser.chainMedian if atomType in analyser.medians.index else 0
        if atomType in analyser.medians.index:
            diffs.append((analyser.medians.loc[atomType]['correctedDensity'] - analyser.chainMedian) / analyser.chainMedian)
            slopes.append(analyser.medians.loc[atomType]['slopes'])
        else:
            diffs.append(0)
            slopes.append(0)
    return diffs, slopes

def main(pdbidfile, resultname, atomType, radius):
    pdbids = []
    with open(pdbidfile, "r") as fileHandleIn:
        for pdbid in fileHandleIn:
            pdbids.append(pdbid[0:4])

    fileHandle = open(resultname, 'w')
    currentAtomType = atomType
    currentRadius = float(radius)
    currentRadii = radiiDefault
    currentRadii[currentAtomType] = currentRadius
    currentSlopes = slopesDefault 
    while True:
        bestMedianDiffs = {}
        bestDiff = bestMedianDiffs[currentAtomType] if bestMedianDiffs else 100 
        while True:
            print("working on", currentAtomType, "with radius", currentRadius, ",", str(datetime.datetime.now()))
            print(currentAtomType, currentRadius, file=fileHandle)
            with multiprocessing.Pool() as pool:
                results = pool.starmap(processFunction, ((pdbid, currentRadii, currentSlopes) for pdbid in pdbids))

            slopes = {}
            diffs = {}
            for atomType in radiiDefault:
                slopes[atomType] = []
                diffs[atomType] = []

            for result in results:
                if result:
                    for i, diff in enumerate(result[0]):
                        if diff: 
                            diffs[sorted(radiiDefault)[i]].append(diff)
                    for i, slope in enumerate(result[1]):
                        if slope != slopesDefault[sorted(slopesDefault)[i]]:
                            slopes[sorted(slopesDefault)[i]].append(slope)

            for key in diffs:
                diffs[key] = np.nanmedian(diffs[key])
            for key in slopes:
                slopes[key] = np.nanmedian(slopes[key])

            print(currentRadii, file=fileHandle)
            print(diffs, file=fileHandle)
            print(slopes, file=fileHandle)
            medianDiff = diffs[currentAtomType]    
            if abs(medianDiff) < abs(bestDiff):
                bestDiff = medianDiff
                bestMedianDiffs = diffs
                bestSlopes = slopes
                bestRadius = currentRadius
                currentRadius = currentRadius + 0.01 if medianDiff < 0 else currentRadius - 0.01
                currentRadii[currentAtomType] = currentRadius
            else: 
                break

        if max(map(abs, bestMedianDiffs.values())) < 0.05:
            break
        else:
            currentAtomType = max(bestMedianDiffs, key=lambda y: abs(bestMedianDiffs[y])) 
            currentRadius = currentRadii[currentAtomType]
            medianDiff = bestMedianDiffs[currentAtomType]
            currentRadius = currentRadius + 0.01 if medianDiff < 0 else currentRadius - 0.01
            currentRadii[currentAtomType] = currentRadius
            currentSlopes = bestSlopes

    fileHandle.close()
    print("Final radii:", currentRadii)


if __name__ == '__main__':
    _, filename, resultname, atomType, radius = sys.argv

    main(filename, resultname, atomType, radius)

