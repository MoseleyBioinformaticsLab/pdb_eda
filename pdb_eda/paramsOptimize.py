#!/usr/bin/python3
"""
pdb_eda radii and slope parameter optimization mode command-line interface
  Optimizes radii and b-factor slopes using a given set of PDB IDs.
  A simple steepest decent optimization approach is utilized.
  This approach is justified by testing and the use of median differences that smooths the error surface.

  Should use a minimum of 1000 PDB entries for optimization.

Usage:
    pdb_eda optimize -h | --help
    pdb_eda optimize <pdbid-file> <log-file> <out-params-file> [options]

Options:
    -h, --help                          Show this screen.
    --atoms=<atoms-file>                Limit optimization to the list of atoms in the JSON atoms-file. This list overrides what is indicated in the params-file. [default: ]
    --sample=<sample-size>              Use a random sample of PDB ids for optimization. [default: 0]
    --max=<max-radius-change>           Maximum to change the radius at each incremental optimization. [default: 0.2]
    --min=<min-radius-change>           Minimum to change the radius at each incremental optimization. [default: 0.005]
    --params=<start-params-file>        Starting parameters filename. [default: ]
    --radius=<start-radius>             Starting radius for the starting atom-type. [default: 0]
    --start=<start-atom-type>           Starting atom type. [default: ]
    --stop=<fractional-difference>      Max fractional difference between atom-specific and chain-specific density conversion allowed for stopping the optimization. [default: 0.05]
    --testing                           Run only a single process for testing purposes.
"""
import os
import sys
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

    paramsFilepath = args["--params"] if args["--params"] else densityAnalysis.paramsPath
    try:
        with open(paramsFilepath, 'r') as jsonFile:
            params = json.load(jsonFile)
            currentRadii = params['radii']
            currentSlopes = params['slopes']
            atoms2Optimize = set(params['optimize']) if 'optimize' in params else []

            if paramsFilepath != densityAnalysis.paramsPath:
                densityAnalysis.setGlobals(params)
    except:
        sys.exit(str("Error: params file \"") + paramsFilename + "\" does not exist or is not parsable.")

    if args["--atoms"]:
        try:
            with open(paramsFilepath, 'r') as jsonFile:
                params = json.load(jsonFile)
                atoms2Optimize = set(params['optimize'])
        except:
            sys.exit(str("Error: atoms file \"") + args["--atoms"] + "\" does not exist or is not parsable.")

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
        (bestMedianDiffs, currentSlopes) = calculateMedianDiffsSlopes(pdbids, params, args["--testing"])
        print("Max Absolute Median Diff: ", max(map(abs, bestMedianDiffs.values())))

        currentAtomType = max(bestMedianDiffs, key=lambda y: abs(bestMedianDiffs[y])) if not args["--start"]  else args["--start"]
        previousRadius = currentRadii[currentAtomType]

        if startingRadius > 0:
            currentRadii[currentAtomType] = startingRadius
        else:
            currentRadii[currentAtomType] = currentRadii[currentAtomType] + radiusIncrement if bestMedianDiffs[currentAtomType] < 0 else currentRadii[currentAtomType] - radiusIncrement

        while True:
            print("Testing", currentAtomType, ": starting radius=", previousRadius, ", new radius=", currentRadii[currentAtomType], ", increment=", radiusIncrement, ", start-time=", str(datetime.datetime.now()))
            print("Testing", currentAtomType, ": starting radius=", previousRadius, ", new radius=", currentRadii[currentAtomType], ", increment=", radiusIncrement, ", start-time=", str(datetime.datetime.now()), file=logFile)

            (medianDiffs, slopes) = calculateMedianDiffsSlopes(pdbids, {**params, "radii" : currentRadii, "slopes" : currentSlopes }, args["--testing"])
            print("Radii: ", currentRadii, file=logFile)
            print("Median Diffs: ", medianDiffs, file=logFile)
            print("Max Absolute Median Diff: ", max(map(abs, medianDiffs.values())))
            print("Max Absolute Median Diff: ", max(map(abs, medianDiffs.values())), file=logFile)
            print("Slopes: ", slopes, file=logFile)

            improved = False
            if abs(medianDiffs[currentAtomType]) < abs(bestMedianDiffs[currentAtomType]):
                bestMedianDiffs = medianDiffs
                currentSlopes = slopes
                improved = True

                try:
                    with open(args["<out-params-file>"] + ".temp", 'w') as jsonFile:
                        print(json.dumps({**params, "radii": currentRadii, "slopes": currentSlopes}, indent=2, sort_keys=True), file=jsonFile)
                except:
                    sys.exit(str("Error: unable to create temporary params file \"") + args["<out-params-file>"] + ".temp" + "\".")
            else:
                currentRadii[currentAtomType] = previousRadius

            testBestMedianDiffs = { atomType:diff for (atomType,diff) in bestMedianDiffs.items() if atomType in atoms2Optimize } if atoms2Optimize else bestMedianDiffs
            maxAtomType = max(testBestMedianDiffs, key=lambda y: abs(testBestMedianDiffs[y]))
            if stoppingFractionalDifference > 0 and max(map(abs, testBestMedianDiffs.values())) < stoppingFractionalDifference:
                break
            elif maxAtomType == currentAtomType:
                if radiusIncrement == minRadiusIncrement and not improved:
                    break

                radiusIncrement = radiusIncrement / 2.0
                if radiusIncrement < minRadiusIncrement:
                    radiusIncrement = minRadiusIncrement

            currentAtomType = maxAtomType
            previousRadius = currentRadii[currentAtomType]
            currentRadii[currentAtomType] = currentRadii[currentAtomType] + radiusIncrement if bestMedianDiffs[currentAtomType] < 0 else currentRadii[currentAtomType] - radiusIncrement

        print("Final Radii: ", currentRadii)
        print("Max Absolute Median Diff", max(map(abs, testBestMedianDiffs.values())))
        outParams = {**params, "radii" : currentRadii, "slopes" : currentSlopes }

    try:
        with open(args["<out-params-file>"], 'w') as jsonFile:
            print(json.dumps(outParams, indent=2, sort_keys=True), file=jsonFile)
    except:
        sys.exit(str("Error: unable to create params file \"") + args["<out-params-file>"] + "\".")

def calculateMedianDiffsSlopes(pdbids, currentParams, testing=False):
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
            results = pool.starmap(processFunction, ((pdbid, currParamsFilename) for pdbid in pdbids))

    diffs = {atomType: [] for atomType in currentParams["radii"]}
    slopes = {atomType: [] for atomType in currentParams["slopes"]}
    for resultFilename in results:
        if resultFilename:
            try:
                with open(resultFilename, 'r') as jsonFile:
                    result = json.load(jsonFile)
                    for atomType, diff in result['diffs'].items():
                        diffs[atomType].append(diff)
                    for atomType, slope in result['slopes'].items():
                        slopes[atomType].append(slope)
                os.remove(resultFilename)
            except:
                pass

    os.remove(currParamsFilename)

    medianDiffs = {key: np.nanmedian(value) for (key, value) in diffs.items()}
    medianSlopes = {key: np.nanmedian(value) for (key, value) in slopes.items()}

    return (medianDiffs, medianSlopes)

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

    analyzer = densityAnalysis.fromPDBid(pdbid)
    if not analyzer:
        return 0 

    analyzer.aggregateCloud(params)
    if not analyzer.densityElectronRatio:
        return 0

    diffs = { atomType:((analyzer.medians['corrected_density_electron_ratio'][atomType] - analyzer.densityElectronRatio) / analyzer.densityElectronRatio) for atomType in params["radii"]
              if atomType in analyzer.medians['corrected_density_electron_ratio'] }
    newSlopes = { atomType:analyzer.medians['slopes'][atomType] for atomType in params["slopes"] if atomType in analyzer.medians['slopes'] }

    resultFilename = createTempJSONFile({ "pdbid" : pdbid, "diffs" : diffs, "slopes" : newSlopes }, "tempResults_")
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

