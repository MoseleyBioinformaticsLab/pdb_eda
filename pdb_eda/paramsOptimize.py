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
    --testing                           Run only a single process for testing purposes.


This mode is often run multiple times using the output parameter file generated in one cycle as the starting parameter file in the next cycle.
Typically, it is good to start with the following series of cycles.
1) --sample=50 --max=0.5 --min=0.01 --stop=0.1 (kill it if it appears to thrash without stopping in a couple hours).
2) --sample=100 --max=0.2 --min=0.005 --stop=0.05  (kill it if appears to thrash after a few hours).
3) --sample=200 --max=0.05 --min=0.001 --stop=0.03  (kill it if appears to thrash after a few hours).
4) --sample=400 --max=0.05 --min=0.001 --stop=0.02  (kill it if appears to thrash after a few hours).
5) use 1000-2000 pdbids --max=0.02 --min=0.001 --stop=0.01  (kill after a day or two it if appears to thrash without stopping).
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
import tempfile
import docopt
import random

from pdb_eda import densityAnalysis

def main():
    args = docopt.docopt(__doc__)
    radiusIncrement = float(args["--max"])
    minRadiusIncrement = float(args["--min"])
    stoppingFractionalDifference = float(args["--stop"])
    startingRadius = float(args["--radius"])
    sampleSize = int(args["--sample"])
    atomsTypes2Optimize = None

    try:
        with open(args["<start-params-file>"], 'r') as jsonFile:
            params = json.load(jsonFile)
            currentRadii = params['radii']
            currentSlopes = params['slopes']
            if not args["--ignore"] and 'optimize' in params:
                atomsTypes2Optimize = set(params['optimize'])

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
        print("Calculating starting median differences: start-time", str(datetime.datetime.now()))
        print("Calculating starting median differences: start-time", str(datetime.datetime.now()), file=logFile)
        (bestMedianDiffs, meanDiffs, overallStdDevDiffs, currentSlopes) = calculateMedianDiffsSlopes(pdbids, params, args["--testing"], args["<pdbid-file>"]+".execution_times")
        print("Max Absolute Median Diff: ", max(map(abs, bestMedianDiffs.values())),
              "Max Abs Diff Mean-Median: ", max([ abs(mean-median) for (mean,median) in zip(bestMedianDiffs.values(), meanDiffs.values())]),
              "Overall Diff StdDev: ", overallStdDevDiffs)
        print("Max Absolute Median Diff: ", max(map(abs, bestMedianDiffs.values())),
              "Max Abs Diff Mean-Median: ", max([ abs(mean-median) for (mean,median) in zip(bestMedianDiffs.values(), meanDiffs.values())]),
              "Overall Diff StdDev: ", overallStdDevDiffs, file=logFile)

        currentAtomType = max(bestMedianDiffs, key=lambda y: abs(bestMedianDiffs[y])) if not args["--start"]  else args["--start"]
        previousRadius = currentRadii[currentAtomType]

        if startingRadius > 0:
            previousDirection = currentRadii[currentAtomType] < startingRadius
            currentRadii[currentAtomType] = startingRadius
        else:
            currentRadii[currentAtomType] = currentRadii[currentAtomType] + radiusIncrement if bestMedianDiffs[currentAtomType] < 0 else currentRadii[currentAtomType] - radiusIncrement
            previousDirection = bestMedianDiffs[currentAtomType] < 0

        while True:
            print("Testing", currentAtomType, ": starting radius=", previousRadius, ", new radius=", currentRadii[currentAtomType], ", increment=", radiusIncrement, ", start-time=", str(datetime.datetime.now()))
            print("Testing", currentAtomType, ": starting radius=", previousRadius, ", new radius=", currentRadii[currentAtomType], ", increment=", radiusIncrement, ", start-time=", str(datetime.datetime.now()), file=logFile)

            (medianDiffs, meanDiffs, overallStdDevDiffs, slopes) = calculateMedianDiffsSlopes(pdbids, {**params, "radii" : currentRadii, "slopes" : currentSlopes }, args["--testing"], args["<pdbid-file>"]+".execution_times")
            print("Radii: ", currentRadii, file=logFile)
            print("Median Diffs: ", medianDiffs, file=logFile)
            print("Max Absolute Median Diff: ", max(map(abs, medianDiffs.values())),
                  "Max Abs Diff Mean-Median: ", max([abs(mean - median) for (mean, median) in zip(medianDiffs.values(), meanDiffs.values())]),
                  "Overall Diff StdDev: ", overallStdDevDiffs)
            print("Max Absolute Median Diff: ", max(map(abs, medianDiffs.values())),
                  "Max Abs Diff Mean-Median: ", max([abs(mean - median) for (mean, median) in zip(medianDiffs.values(), meanDiffs.values())]),
                  "Overall Diff StdDev: ", overallStdDevDiffs, file=logFile)
            print("Slopes: ", slopes, file=logFile)

            improved = False
            if abs(medianDiffs[currentAtomType]) < abs(bestMedianDiffs[currentAtomType]):
                bestMedianDiffs = medianDiffs
                currentSlopes = slopes
                improved = True

                print("       ", currentAtomType, "New Radius Accepted: ", currentRadii[currentAtomType])
                print("       ", currentAtomType, "New Radius Accepted: ", currentRadii[currentAtomType], file=logFile)

                try:
                    with open(args["<out-params-file>"] + ".temp", 'w') as jsonFile:
                        print(json.dumps({**params, "radii": currentRadii, "slopes": currentSlopes}, indent=2, sort_keys=True), file=jsonFile)
                except:
                    sys.exit(str("Error: unable to create temporary params file \"") + args["<out-params-file>"] + ".temp" + "\".")
            else:
                print("       ", currentAtomType, "New Radius Rejected: ", currentRadii[currentAtomType])
                print("       ", currentAtomType, "New Radius Rejected: ", currentRadii[currentAtomType], file=logFile)
                currentRadii[currentAtomType] = previousRadius

            testBestMedianDiffs = { atomType:diff for (atomType,diff) in bestMedianDiffs.items() if atomType in atomsTypes2Optimize } if atomsTypes2Optimize else bestMedianDiffs
            maxAtomType = max(testBestMedianDiffs, key=lambda y: abs(testBestMedianDiffs[y]))
            if stoppingFractionalDifference > 0 and max(map(abs, testBestMedianDiffs.values())) < stoppingFractionalDifference:
                break
            elif maxAtomType == currentAtomType and (not improved or previousDirection != (bestMedianDiffs[currentAtomType] < 0)):
                if radiusIncrement == minRadiusIncrement:
                    break

                radiusIncrement = radiusIncrement / 2.0
                if radiusIncrement < minRadiusIncrement:
                    radiusIncrement = minRadiusIncrement

            currentAtomType = maxAtomType
            previousRadius = currentRadii[currentAtomType]
            currentRadii[currentAtomType] = currentRadii[currentAtomType] + radiusIncrement if bestMedianDiffs[currentAtomType] < 0 else currentRadii[currentAtomType] - radiusIncrement
            previousDirection = bestMedianDiffs[currentAtomType] < 0

        print("Final Radii: ", currentRadii)
        print("Max Absolute Median Diff", max(map(abs, testBestMedianDiffs.values())))
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
    currParamsFilename = createTempJSONFile(currentParams, "tempParams_")

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

    medianDiffs = {key: np.nanmedian(value) for (key, value) in diffs.items()}
    meanDiffs = {key: np.mean(value) for (key, value) in diffs.items()}
    overallStdDevDiffs = np.std([item for values in diffs.values() for item in values])
    medianSlopes = {key: np.nanmedian(value) for (key, value) in slopes.items()}

    return (medianDiffs, meanDiffs, overallStdDevDiffs, medianSlopes)

def processFunction(pdbid, paramsFilepath):
    """Process function to analyze a single pdb entry.

    :param :py:class:`str` pdbid: pdbid for entry to download and analyze.
    :param :py:class:`str` paramsFilepath: filepath to the radii and slopes parameters.
    :return: resultFilename
    :rtype: :py:class:`str` or 0
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
              if atomType in analyzer.medians['corrected_density_electron_ratio'] }
    newSlopes = { atomType:analyzer.medians['slopes'][atomType] for atomType in params["slopes"] if atomType in analyzer.medians['slopes'] }

    elapsedTime = time.process_time() - startTime
    resultFilename = createTempJSONFile({ "pdbid" : pdbid, "diffs" : diffs, "slopes" : newSlopes, 'resolution' : analyzer.pdbObj.header.resolution, "execution_time" : elapsedTime }, "tempResults_")
    analyzer=0
    gc.collect()
    return resultFilename

def createTempJSONFile(data, filenamePrefix):
    """Creates a temporary JSON file and returns its filename.

    :param data:  data to save into the JSON file.
    :param :py:class:`str` filenamePrefix: temporary filename prefix.
    :return: filename
    :rtype: :py:class:`str`
    """
    dirname = os.getcwd()
    filename = 0
    with tempfile.NamedTemporaryFile(mode='w', buffering=1, dir=dirname, prefix=filenamePrefix, delete=False) as tempFile:
        json.dump(data,tempFile)
        filename = tempFile.name
    return filename

