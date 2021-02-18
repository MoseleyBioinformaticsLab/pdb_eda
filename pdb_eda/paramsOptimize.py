#!/usr/bin/python3
"""
pdb_eda radii and slope parameter optimization mode command-line interface
  Optimizes radii and b-factor slopes using a given set of PDB IDs.
  A simple steepest decent optimization approach is utilized.
  This approach is justified by testing and the use of median differences that smooths the error surface.

  Should use a minimum of 1000 PDB entries for optimization.

Usage:
    pdb_eda optimize -h | --help
    pdb_eda optimize <start-params-file> <pdbid-file> <log-file> <out-params-file> [options]

Options:
    -h, --help                          Show this screen.
    --ignore                            Ignore the "optimize" atom type limit in the parameter file.
    --sample=<sample-size>              Use a random sample of PDB ids for optimization. [default: 0]
    --max=<max-radius-change>           Maximum to change the radius at each incremental optimization. [default: 0.2]
    --min=<min-radius-change>           Minimum to change the radius at each incremental optimization. [default: 0.001]
    --radius=<start-radius>             Starting radius for the starting atom-type. [default: 0]
    --start=<start-atom-type>           Starting atom type. [default: ]
    --stop=<fractional-difference>      Max fractional difference between atom-specific and chain-specific density conversion allowed for stopping the optimization. [default: 0.02]
    --unweighted                        Pick atom type to optimize next without weighting based on occurrence across PDB entries.
    --testing                           Run only a single process for testing purposes.


This mode is often run multiple times using the output parameter file generated in one cycle as the starting parameter file in the next cycle.
Typically, it is good to start with the following series of cycles.
1) --sample=50 --max=0.2 --min=0.01 --stop=0.1 (kill it if it appears to thrash without stopping in a couple hours).
2) --sample=100 --max=0.2 --min=0.005 --stop=0.05  (kill it if appears to thrash after a few hours).
3) --sample=200 --max=0.05 --min=0.001 --stop=0.03  (kill it if appears to thrash after a few hours).
4) --sample=400 --max=0.05 --min=0.001 --stop=0.02  (kill it if appears to thrash after a few hours).
5) (use 1000-2000 pdbids) --max=0.02 --min=0.001 --stop=0.015  (kill after a day or two it if appears to thrash without stopping).
If you need to kill the run, don't worry, there is a temp output parameter file with the last improvement.
"""
import os
import gc
import sys
import time
import json
import numpy as np
import multiprocessing
import datetime
import docopt
import random

from . import densityAnalysis
from . import fileUtils
from . import __version__

def main():
    args = docopt.docopt(__doc__, version=__version__)
    maxRadiusIncrement = float(args["--max"])
    radiusIncrement = maxRadiusIncrement
    minRadiusIncrement = float(args["--min"])
    stoppingFractionalDifference = float(args["--stop"])
    startingRadius = float(args["--radius"])
    sampleSize = int(args["--sample"])
    atomTypes2Optimize = None

    try:
        with open(args["<start-params-file>"], 'r') as jsonFile:
            params = json.load(jsonFile)
            currentRadii = params['radii']
            currentSlopes = params['slopes']
            if not args["--ignore"] and 'optimize' in params:
                atomTypes2Optimize = set(params['optimize'])

            densityAnalysis.setGlobals(params)
    except:
        sys.exit(str("Error: params file \"") + args["<start-params-file>"] + "\" does not exist or is not parsable.")

    if args["--start"] != "" and args["--start"] not in currentRadii:
        sys.exit(str("Error: starting atom \"") + args["--start"] + "\" is not valid.")

    try:
        pdbids = []
        with open(args["<pdbid-file>"], "r") as textFile:
            for pdbid in textFile:
                pdbids.append(pdbid[0:4])
    except:
        sys.exit(str("Error: PDB IDs file \"") + args["<pdbid-file>"] + "\" does not exist or is not parsable.")

    if sampleSize > 0:
        pdbids = random.sample(pdbids,sampleSize)

    with open(args["<log-file>"], 'w') as logFile:
        print(args, file=logFile)

        print("PDB IDs:",",".join(pdbids), file=logFile)
        print("Calculating start median differences: start-time", str(datetime.datetime.now()))
        print("Calculating start median differences: start-time", str(datetime.datetime.now()), file=logFile)

        (bestMedianDiffs, meanDiffs, overallStdDevDiffs, currentSlopes, sizes) = calculateMedianDiffsSlopes(pdbids, params, args["--testing"], args["<pdbid-file>"]+".execution_times")

        maxSize = max([sizes[atomType] for atomType in bestMedianDiffs.keys() if not atomTypes2Optimize or atomType in atomTypes2Optimize])
        print("Radii:", currentRadii, file=logFile)
        print("Median Diffs:", bestMedianDiffs, file=logFile)
        print("Max Absolute Weighted Median Diff:", max([abs(bestMedianDiffs[atomType] * sizes[atomType] / maxSize) for atomType in bestMedianDiffs.keys() if not atomTypes2Optimize or atomType in atomTypes2Optimize]),
              "Weighted Diff StdDev:", overallStdDevDiffs,
              "Max Size:", maxSize)
        print("Max Absolute Weighted Median Diff:", max([abs(bestMedianDiffs[atomType] * sizes[atomType] / maxSize) for atomType in bestMedianDiffs.keys() if not atomTypes2Optimize or atomType in atomTypes2Optimize]),
              "Weighted Diff StdDev:", overallStdDevDiffs,
              "Max Size:", maxSize, file=logFile)
        print("Max Absolute Median Diff:", max([abs(bestMedianDiffs[atomType]) for atomType in bestMedianDiffs.keys() if not atomTypes2Optimize or atomType in atomTypes2Optimize]),
              "Max Abs Diff Mean-Median:", max([abs(meanDiffs[atomType] - bestMedianDiffs[atomType]) for atomType in bestMedianDiffs.keys() if not atomTypes2Optimize or atomType in atomTypes2Optimize]),
              "Mean Abs Diff Mean-Median:", np.mean([abs(meanDiffs[atomType] - bestMedianDiffs[atomType]) for atomType in bestMedianDiffs.keys() if not atomTypes2Optimize or atomType in atomTypes2Optimize]))
        print("Max Absolute Median Diff:", max([abs(bestMedianDiffs[atomType]) for atomType in bestMedianDiffs.keys() if not atomTypes2Optimize or atomType in atomTypes2Optimize]),
              "Max Abs Diff Mean-Median:", max([abs(meanDiffs[atomType] - bestMedianDiffs[atomType]) for atomType in bestMedianDiffs.keys() if not atomTypes2Optimize or atomType in atomTypes2Optimize]),
              "Mean Abs Diff Mean-Median:", np.mean([abs(meanDiffs[atomType] - bestMedianDiffs[atomType]) for atomType in bestMedianDiffs.keys() if not atomTypes2Optimize or atomType in atomTypes2Optimize]), file=logFile)

        testBestMedianDiffs = {atomType: diff for (atomType, diff) in bestMedianDiffs.items() if atomType in atomTypes2Optimize} if atomTypes2Optimize else bestMedianDiffs
        if args["--unweighted"]:
            currentAtomType = max(testBestMedianDiffs, key=lambda y: abs(testBestMedianDiffs[y])) if not args["--start"]  else args["--start"]
        else:
            currentAtomType = max(testBestMedianDiffs, key=lambda y: abs(testBestMedianDiffs[y]) * sizes[y]) if not args["--start"] else args["--start"]
        previousRadius = currentRadii[currentAtomType]

        if startingRadius > 0:
            previousDirection = currentRadii[currentAtomType] < startingRadius
            currentRadii[currentAtomType] = startingRadius
        else:
            currentRadii[currentAtomType] = currentRadii[currentAtomType] + radiusIncrement if bestMedianDiffs[currentAtomType] < 0 else currentRadii[currentAtomType] - radiusIncrement
            previousDirection = bestMedianDiffs[currentAtomType] < 0

        estimatedRadiusIncrement = {atomType:0 for atomType in currentRadii.keys()}
        while True:
            print("Testing ", currentAtomType, ": starting radius=", previousRadius, ", new radius=", currentRadii[currentAtomType],
                  ", current weighted median difference=",bestMedianDiffs[currentAtomType] * sizes[currentAtomType] / maxSize, str("(")+str(bestMedianDiffs[currentAtomType])+str(")"), ", size=",sizes[currentAtomType])
            print("Testing ", currentAtomType, ": starting radius=", previousRadius, ", new radius=", currentRadii[currentAtomType],
                  ", current median difference=",bestMedianDiffs[currentAtomType], str("(")+str(bestMedianDiffs[currentAtomType])+str(")"), ", size=",sizes[currentAtomType], file=logFile)
            print("Calculating next  median differences: start-time", str(datetime.datetime.now()))
            print("Calculating next  median differences: start-time", str(datetime.datetime.now()), file=logFile)

            (medianDiffs, meanDiffs, overallStdDevDiffs, slopes, sizes) = calculateMedianDiffsSlopes(pdbids, {**params, "radii" : currentRadii, "slopes" : currentSlopes }, args["--testing"],
                                                                                                     args["<pdbid-file>"]+".execution_times")

            maxSize = max([sizes[atomType] for atomType in medianDiffs.keys() if not atomTypes2Optimize or atomType in atomTypes2Optimize])
            print("Radii:", currentRadii, file=logFile)
            print("Median Diffs:", medianDiffs, file=logFile)
            print("Max Absolute Weighted Median Diff:", max([abs(medianDiffs[atomType] * sizes[atomType] / maxSize) for atomType in medianDiffs.keys() if not atomTypes2Optimize or atomType in atomTypes2Optimize]),
                  "Weighted Diff StdDev:", overallStdDevDiffs,
                  "Max Size:", maxSize)
            print("Max Absolute Weighted Median Diff:", max([abs(medianDiffs[atomType] * sizes[atomType] / maxSize) for atomType in medianDiffs.keys() if not atomTypes2Optimize or atomType in atomTypes2Optimize]),
                  "Weighted Diff StdDev:", overallStdDevDiffs,
                  "Max Size:", maxSize, file=logFile)
            print("Max Absolute Median Diff:", max([abs(medianDiffs[atomType]) for atomType in medianDiffs.keys() if not atomTypes2Optimize or atomType in atomTypes2Optimize]),
                  "Max Abs Diff Mean-Median:", max([abs(meanDiffs[atomType] - medianDiffs[atomType]) for atomType in medianDiffs.keys() if not atomTypes2Optimize or atomType in atomTypes2Optimize]),
                  "Mean Abs Diff Mean-Median:", np.mean([abs(meanDiffs[atomType] - medianDiffs[atomType]) for atomType in medianDiffs.keys() if not atomTypes2Optimize or atomType in atomTypes2Optimize]))
            print("Max Absolute Median Diff:", max([abs(medianDiffs[atomType]) for atomType in medianDiffs.keys() if not atomTypes2Optimize or atomType in atomTypes2Optimize]),
                  "Max Abs Diff Mean-Median:", max([abs(meanDiffs[atomType] - medianDiffs[atomType]) for atomType in medianDiffs.keys() if not atomTypes2Optimize or atomType in atomTypes2Optimize]),
                  "Mean Abs Diff Mean-Median:", np.mean([abs(meanDiffs[atomType] - medianDiffs[atomType]) for atomType in medianDiffs.keys() if not atomTypes2Optimize or atomType in atomTypes2Optimize]), file=logFile)
            print("Slopes:", slopes, file=logFile)

            improved = False
            if abs(medianDiffs[currentAtomType]) <= abs(bestMedianDiffs[currentAtomType]):
                if abs(medianDiffs[currentAtomType]) < abs(bestMedianDiffs[currentAtomType]) and previousDirection == (medianDiffs[currentAtomType] < 0):
                    estimatedRadiusIncrement[currentAtomType] = 0.80 * (currentRadii[currentAtomType] - previousRadius) * medianDiffs[currentAtomType] / (bestMedianDiffs[currentAtomType] - medianDiffs[currentAtomType])
                else:
                    estimatedRadiusIncrement[currentAtomType] = 0
                bestMedianDiffs = medianDiffs
                currentSlopes = slopes
                improved = True if abs(medianDiffs[currentAtomType]) < abs(bestMedianDiffs[currentAtomType]) else 2

                print("Accepted", currentAtomType, ": new radius=", currentRadii[currentAtomType],", current weighted median difference=",bestMedianDiffs[currentAtomType] * sizes[currentAtomType] / maxSize,
                      str("(")+str(bestMedianDiffs[currentAtomType])+str(")"), ", size=",sizes[currentAtomType])
                print("Accepted", currentAtomType, ": new radius=", currentRadii[currentAtomType],", current weighted median difference=",bestMedianDiffs[currentAtomType] * sizes[currentAtomType] / maxSize,
                      str("(")+str(bestMedianDiffs[currentAtomType])+str(")"), ", size=",sizes[currentAtomType], file=logFile)

                try:
                    with open(args["<out-params-file>"] + ".temp", 'w') as jsonFile:
                        print(json.dumps({**params, "radii": currentRadii, "slopes": currentSlopes}, indent=2, sort_keys=True), file=jsonFile)
                except:
                    sys.exit(str("Error: unable to create temporary params file \"") + args["<out-params-file>"] + ".temp" + "\".")
            else:
                estimatedRadiusIncrement[currentAtomType] = 0
                print("Rejected", currentAtomType, ": new radius=", currentRadii[currentAtomType])
                print("Rejected", currentAtomType, ": new radius=", currentRadii[currentAtomType], file=logFile)
                currentRadii[currentAtomType] = previousRadius

            testBestMedianDiffs = { atomType:diff for (atomType,diff) in bestMedianDiffs.items() if atomType in atomTypes2Optimize } if atomTypes2Optimize else bestMedianDiffs
            if args["--unweighted"]:
                maxAtomType = max(testBestMedianDiffs, key=lambda y: abs(testBestMedianDiffs[y]))
            else:
                maxAtomType = max(testBestMedianDiffs, key=lambda y: abs(testBestMedianDiffs[y]) * sizes[y])
            if stoppingFractionalDifference > 0 and max([abs(value * sizes[atomType] / maxSize) for atomType,value in testBestMedianDiffs.items()]) < stoppingFractionalDifference:
                break
            elif maxAtomType == currentAtomType:
                if not improved or previousDirection != (bestMedianDiffs[currentAtomType] < 0):
                    if radiusIncrement == minRadiusIncrement:
                        break

                    radiusIncrement = radiusIncrement / 2.0
                    if radiusIncrement < minRadiusIncrement:
                        radiusIncrement = minRadiusIncrement
                elif improved == 2:
                    radiusIncrement = radiusIncrement * 1.5
                    if radiusIncrement > maxRadiusIncrement:
                        radiusIncrement = maxRadiusIncrement

                print("New Radius Increment:", radiusIncrement)
                print("New Radius Increment:", radiusIncrement, file=logFile)

            currentAtomType = maxAtomType
            previousRadius = currentRadii[currentAtomType]
            if abs(estimatedRadiusIncrement[currentAtomType]) > 0:
                currentRadii[currentAtomType] = currentRadii[currentAtomType] + estimatedRadiusIncrement[currentAtomType]
            else:
                currentRadii[currentAtomType] = currentRadii[currentAtomType] + radiusIncrement if bestMedianDiffs[currentAtomType] < 0 else currentRadii[currentAtomType] - radiusIncrement
            previousDirection = bestMedianDiffs[currentAtomType] < 0
            gc.collect() # force garbage collection to decrease memory use.

        print("Final Radii:", currentRadii)
        print("Max Absolute Weighted Median Diff:", max([abs(value * sizes[atomType] / maxSize) for atomType,value in testBestMedianDiffs.items()]))
        outParams = {**params, "radii" : currentRadii, "slopes" : currentSlopes }

    try:
        with open(args["<out-params-file>"], 'w') as jsonFile:
            print(json.dumps(outParams, indent=2, sort_keys=True), file=jsonFile)
    except:
        sys.exit(str("Error: unable to create params file \"") + args["<out-params-file>"] + "\".")

def calculateMedianDiffsSlopes(pdbids, currentParams, testing=False, executionTimesFilename=None):
    """Calculates the median diffs and slopes across a list of pdb entries.

    :param :py:class:`list` pdbids: list of pdbids to process.
    :param :py:class:`dict` currentParams:  parameters.
    :return: diffs_slopes_tuple
    :rtype: :py:class:`tuple`
    """
    currParamsFilename = fileUtils.createTempJSONFile(currentParams, "tempParams_")

    if testing:
        results = [processFunction(pdbid) for pdbid in pdbids]
    else:
        with multiprocessing.Pool() as pool:
            results = pool.starmap(processFunction, ((pdbid, currParamsFilename) for pdbid in pdbids), chunksize=1)

    diffs = {atomType: [] for atomType in currentParams["radii"]}
    slopes = {atomType: [] for atomType in currentParams["slopes"]}
    executionTimes = {}
    for resultFilename in results:
        if resultFilename:
            try:
                with open(resultFilename, 'r') as jsonFile:
                    result = json.load(jsonFile)
                    for atomType, diff in result['diffs'].items():
                        diffs[atomType].append(diff)
                    for atomType, slope in result['slopes'].items():
                        slopes[atomType].append(slope)
                    executionTimes[result["pdbid"]] = result["execution_time"]
                os.remove(resultFilename)
            except:
                pass

    os.remove(currParamsFilename)

    # put longest running jobs at the beginning to improve the running time of each iteration.
    pdbids.sort(key=lambda x : executionTimes[x] if x in executionTimes else 0, reverse=True)

    # print execution times
    if executionTimesFilename:
        with open(executionTimesFilename, "w") as txtFile:
            print("\n".join(pdbid + "  - " + str(executionTimes[pdbid] if pdbid in executionTimes else 0) for pdbid in pdbids), file=txtFile)

    medianDiffs = {key: (np.nanmedian(value) if (value and not np.isnan(value).all()) else 0) for (key, value) in diffs.items() }
    meanDiffs = {key: (np.nanmean(value) if (value and not np.isnan(value).all()) else 0) for (key, value) in diffs.items() }
    sizeDiffs = {key: sum(~np.isnan(value)) for (key, value) in diffs.items() }
    squaredDiffs = [item ** 2 for values in diffs.values() for item in values if not np.isnan(item)]
    overallStdDevDiffs = np.sqrt(sum(squaredDiffs)/(len(squaredDiffs)-1))
    medianSlopes = {key: np.nanmedian(value) for (key, value) in slopes.items()}

    return (medianDiffs, meanDiffs, overallStdDevDiffs, medianSlopes, sizeDiffs)

def processFunction(pdbid, paramsFilepath):
    """Process function to analyze a single pdb entry.

    :param :py:class:`str` pdbid: pdbid for entry to download and analyze.
    :param :py:class:`str` paramsFilepath: filepath to the radii and slopes parameters.
    :return: resultFilename or 0
    :rtype: :py:class:`str` or :py:class:`int`
    """
    try:
        with open(paramsFilepath, 'r') as jsonFile:
            params = json.load(jsonFile)
    except:
        return 0

    startTime = time.process_time()

    analyzer = densityAnalysis.fromPDBid(pdbid)
    if not analyzer:
        return 0 

    analyzer.aggregateCloud(params)
    if not analyzer.densityElectronRatio:
        return 0

    diffs = { atomType:((analyzer.medians['corrected_density_electron_ratio'][atomType] - analyzer.densityElectronRatio) / analyzer.densityElectronRatio) for atomType in params["radii"]
              if atomType in analyzer.medians['corrected_density_electron_ratio'] and not np.isnan(analyzer.medians['corrected_density_electron_ratio'][atomType]) }
    newSlopes = { atomType:analyzer.medians['slopes'][atomType] for atomType in params["slopes"] if atomType in analyzer.medians['slopes'] and not np.isnan(analyzer.medians['slopes'][atomType]) }

    elapsedTime = time.process_time() - startTime
    resultFilename = fileUtils.createTempJSONFile({ "pdbid" : pdbid, "diffs" : diffs, "slopes" : newSlopes, 'resolution' : analyzer.pdbObj.header.resolution, "execution_time" : elapsedTime }, "tempResults_")

    # force garbage collection to decrease memory use.
    analyzer=0
    diffs=0
    newSlopes=0
    gc.collect()

    return resultFilename


