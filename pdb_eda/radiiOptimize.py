#!/usr/bin/python3
"""
optimizeParams.py
  Optimizes radii and b-factor slopes using a given set of PDB IDs.

Usage:
    optimizeParams -h | --help
    optimizeParams <pdbid-file> <log-file> <final-params-file> <start-atom-type> [options]

Options:
    -h, --help                          Show this screen.
    --params=<start-params-file>        Starting params file. [default: ""]
    --radius=<start-radius>             Starting radius for the starting atom-type. [default: 0]
    --change=<radius-change>            How much to change the radius at each incremental optimization. [default: 0.01]
    --stop=<fractional-difference>      Max fractional difference between atom-specific and chain-specific density conversion allowed for stopping the optimization. [default: 0.05]
"""
import os
import sys
import json
import numpy as np
import multiprocessing
import datetime
import tempfile
import docopt

from pdb_eda import densityAnalysis

defaultParamsFilename = os.path.join(os.path.dirname(__file__), 'conf/intermediate_radii_slope_param.json')

def processFunction(pdbid, paramsPath):
    try:
        with open(paramsPath, 'r') as jsonFile:
            params = json.load(jsonFile)
            radii = params['radii']
            slopes = params['slopes']
    except:
        return 0

    analyser = densityAnalysis.fromPDBid(pdbid)
    if not analyser:
        return 0 

    analyser.aggregateCloud(radii, slopes)
    if not analyser.chainMedian:
        return 0

    diffs = []
    newSlopes = []
    for atomType in sorted(radii):
        if atomType in analyser.medians.index:
            diffs.append((analyser.medians.loc[atomType]['correctedDensity'] - analyser.chainMedian) / analyser.chainMedian)
            newSlopes.append(analyser.medians.loc[atomType]['slopes'])
        else:
            diffs.append(0)
            newSlopes.append(0)

    resultFilename = createTempJSONFile({ "pdbid" : pdbid, "diffs" : diffs, "slopes" : newSlopes }, "tempResults_")
    return resultFilename

def createTempJSONFile(data, filenamePrefix):
    dirname = os.getcwd()
    filename = 0
    with tempfile.NamedTemporaryFile(mode='w', buffering=1, dir=dirname, prefix=filenamePrefix, delete=False) as tempFile:
        json.dump(data,tempFile)
        filename = tempFile.name
    return filename


def main(args):
    radiusIncrement = float(args["--change"])
    stoppingFractionalDifference = float(args["--stop"])
    startingRadius = float(args["--radius"])

    paramsFilename = args["--params"] if args["--params"] else defaultParamsFilename
    try:
        with open(paramsFilename, 'r') as jsonFile:
            params = json.load(jsonFile)
            radiiDefault = params['radii']
            slopesDefault = params['slopes']
    except:
        sys.exit(str("Error: params file \"") + paramsFilename + "\" does not exist or is not parsable.")

    try:
        pdbids = []
        with open(args["<pdbid-file>"], "r") as textFile:
            for pdbid in textFile:
                pdbids.append(pdbid[0:4])
    except:
        sys.exit(str("Error: PDB IDs file \"") + args["<pdbid-file>"] + "\" does not exist or is not parsable.")

    try:
        with open(args["<log-file>"], 'w') as logFile:
            print(args, file=logFile)
            currentAtomType = args["<start-atom-type>"]
            currentRadii = radiiDefault
            if startingRadius > 0:
                currentRadii[currentAtomType] = startingRadius
            currentRadius = currentRadii[currentAtomType]
            previousRadius = currentRadius
            currentSlopes = slopesDefault
            bestSlopes = currentSlopes
            bestDiff = 100
            while True:
                while True:
                    print("working on", currentAtomType, "with radius", currentRadius, ",", str(datetime.datetime.now()))
                    print(currentAtomType, currentRadius, file=logFile)
                    currParamsFilename = createTempJSONFile({"radii": currentRadii, "slopes": currentSlopes}, "tempParams_")

                    with multiprocessing.Pool() as pool:
                        results = pool.starmap(processFunction, ((pdbid, currParamsFilename) for pdbid in pdbids))

                    slopes = { atomType:[] for atomType in radiiDefault }
                    diffs = { atomType:[] for atomType in radiiDefault }
                    for resultFilename in results:
                        if resultFilename:
                            try:
                                with open(resultFilename, 'r') as jsonFile:
                                    result = json.load(jsonFile)
                                    for i, diff in enumerate(result['diffs']):
                                        if diff:
                                            diffs[sorted(radiiDefault)[i]].append(diff)
                                    for i, slope in enumerate(result['slopes']):
                                        if slope != slopesDefault[sorted(slopesDefault)[i]]:
                                            slopes[sorted(slopesDefault)[i]].append(slope)
                                os.remove(resultFilename)
                            except:
                                pass

                    diffs = { key:np.nanmedian(value) for (key,value) in diffs.items() }
                    slopes = { key:np.nanmedian(value) for (key,value) in slopes.items() }

                    print(currentRadii, file=logFile)
                    print(diffs, file=logFile)
                    print(slopes, file=logFile)
                    medianDiff = diffs[currentAtomType]
                    if abs(medianDiff) < abs(bestDiff):
                        bestDiff = medianDiff
                        bestMedianDiffs = diffs
                        bestSlopes = slopes
                        currentSlopes = slopes # may be dangerous, because slopes are changing one step behind currentRadius, but the alternative is worse.
                        previousRadius = currentRadius
                        currentRadius = currentRadius + radiusIncrement if medianDiff < 0 else currentRadius - radiusIncrement
                        currentRadii[currentAtomType] = currentRadius
                    else:
                        currentRadii[currentAtomType] = previousRadius
                        break

                os.remove(currParamsFilename) # Remove unneeded params file.

                maxAtomType = max(bestMedianDiffs, key=lambda y: abs(bestMedianDiffs[y]))
                if max(map(abs, bestMedianDiffs.values())) < stoppingFractionalDifference or maxAtomType == currentAtomType:
                    break
                else:
                    currentAtomType = max(bestMedianDiffs, key=lambda y: abs(bestMedianDiffs[y]))
                    currentRadius = currentRadii[currentAtomType]
                    previousRadius = currentRadius
                    medianDiff = bestMedianDiffs[currentAtomType]
                    bestDiff = medianDiff
                    currentRadius = currentRadius + radiusIncrement if medianDiff < 0 else currentRadius - radiusIncrement
                    currentRadii[currentAtomType] = currentRadius
                    currentSlopes = bestSlopes # likely redundant now.
    except:
        sys.exit(str("Error: unable to open log file \"") + args["<log-file>"] + "\".")

    print("Final radii:", currentRadii)

    try:
        with open(args["<final-params-file>"], 'w') as jsonFile:
            json.dump({ "radii" : currentRadii, "slopes" : currentSlopes },jsonFile)
    except:
        sys.exit(str("Error: unable to create final params file \"") + args["<final-params-file>"] + "\".")



if __name__ == '__main__':
    args = docopt.docopt(__doc__)
    main(args)

