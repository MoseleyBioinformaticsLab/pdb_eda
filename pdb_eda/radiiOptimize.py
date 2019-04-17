import os
import sys
import json
import numpy as np
import multiprocessing
import datetime

from pdb_eda import densityAnalysis

## Radii and slopes from an intial analysis on 100 structures
radiiParamPath = os.path.join(os.path.dirname(__file__), 'conf/intermediate_radii_slope_param.json')
with open(radiiParamPath, 'r') as fh:
    radiiParams = json.load(fh)

radiiDefault = radiiParams['radii']
slopesDefault = radiiParams['slopes']

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

