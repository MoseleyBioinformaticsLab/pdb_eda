#!/usr/bin/python3
"""
pdb_eda radii and slope parameter optimization mode command-line interface
  Optimizes radii and b-factor slopes using a given set of PDB IDs.

Usage:
    pdb_eda optimize -h | --help
    pdb_eda optimize <pdbid-file> <log-file> <final-params-file> [options]

Options:
    -h, --help                          Show this screen.
    --params=<start-params-file>        Starting params file. [default: ]
    --radius=<start-radius>             Starting radius for the starting atom-type. [default: 0]
    --atom=<start-atom-type>            Starting atom type. [default: ]
    --max=<max-radius-change>           Maximum to change the radius at each incremental optimization. [default: 0.5]
    --min=<min-radius-change>           Minimum to change the radius at each incremental optimization. [default: 0.005]
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

defaultParamsFilepath = os.path.join(os.path.dirname(__file__), 'conf/intermediate_radii_slope_param.json')

def main():
    args = docopt.docopt(__doc__)
    radiusIncrement = float(args["--max"])
    minRadiusIncrement = float(args["--min"])
    stoppingFractionalDifference = float(args["--stop"])
    startingRadius = float(args["--radius"])

    paramsFilepath = args["--params"] if args["--params"] else defaultParamsFilepath
    try:
        with open(paramsFilepath, 'r') as jsonFile:
            params = json.load(jsonFile)
            currentRadii = params['radii']
            currentSlopes = params['slopes']
    except:
        sys.exit(str("Error: params file \"") + paramsFilename + "\" does not exist or is not parsable.")

    if args["--atom"] != "" and args["--atom"] not in currentRadii:
        sys.exit(str("Error: starting atom \"") + args["--atom"] + "\" is not valid.")

    try:
        pdbids = []
        with open(args["<pdbid-file>"], "r") as textFile:
            for pdbid in textFile:
                pdbids.append(pdbid[0:4])
    except:
        sys.exit(str("Error: PDB IDs file \"") + args["<pdbid-file>"] + "\" does not exist or is not parsable.")

    with open(args["<log-file>"], 'w') as logFile:
        print(args, file=logFile)

        print("Calculating starting median differences: ", str(datetime.datetime.now()))
        (bestMedianDiffs, currentSlopes) = calculateMedianDiffsSlopes(pdbids, currentRadii, currentSlopes)

        currentAtomType = max(bestMedianDiffs, key=lambda y: abs(bestMedianDiffs[y])) if not args["--atom"]  else args["--atom"]
        previousRadius = currentRadii[currentAtomType]

        if startingRadius > 0:
            currentRadii[currentAtomType] = startingRadius
        else:
            currentRadii[currentAtomType] = currentRadii[currentAtomType] + radiusIncrement if bestMedianDiffs[currentAtomType] < 0 else currentRadii[currentAtomType] - radiusIncrement

        while True:
            print("Working on ", currentAtomType, "with radius ", currentRadii[currentAtomType], " with increment ", radiusIncrement, ": ", str(datetime.datetime.now()))
            print(currentAtomType, currentRadii[currentAtomType], file=logFile)

            (medianDiffs, slopes) = calculateMedianDiffsSlopes(pdbids, currentRadii, currentSlopes)
            print("Radii: ", currentRadii, file=logFile)
            print("Median Diffs: ", medianDiffs, file=logFile)
            print("Max Absolute Median Diff: ", max(map(abs, medianDiffs.values())))
            print("Slopes: ", slopes, file=logFile)

            if abs(medianDiffs[currentAtomType]) < abs(bestMedianDiffs[currentAtomType]):
                bestMedianDiffs = medianDiffs
                currentSlopes = slopes
            else:
                currentRadii[currentAtomType] = previousRadius

            maxAtomType = max(bestMedianDiffs, key=lambda y: abs(bestMedianDiffs[y]))
            if stoppingFractionalDifference > 0 and max(map(abs, bestMedianDiffs.values())) < stoppingFractionalDifference:
                break
            elif maxAtomType == currentAtomType:
                radiusIncrement = radiusIncrement / 2.0

                if radiusIncrement < minRadiusIncrement:
                    break

            currentAtomType = maxAtomType
            previousRadius = currentRadii[currentAtomType]
            currentRadii[currentAtomType] = currentRadii[currentAtomType] + radiusIncrement if bestMedianDiffs[currentAtomType] < 0 else currentRadii[currentAtomType] - radiusIncrement

    print("Final Radii: ", currentRadii)
    print("Max Absolute Median Diff", max(map(abs, bestMedianDiffs.values())))

    try:
        with open(args["<final-params-file>"], 'w') as jsonFile:
            json.dump({ "radii" : currentRadii, "slopes" : currentSlopes },jsonFile)
    except:
        sys.exit(str("Error: unable to create final params file \"") + args["<final-params-file>"] + "\".")

def calculateMedianDiffsSlopes(pdbids, currentRadii, currentSlopes):
    """Calculates the median diffs and slopes across a list of pdb entries.

    :param :py:class:`list` pdbids: list of pdbids to process.
    :param :py:class:`dict` currentRadii:  atom type radii to use.
    :param :py:class:`dict` currentSlopes: atom type b-factor slopes to use.
    :return: diffs_slopes_tuple
    :rtype: :py:class:`tuple`
    """
    currParamsFilename = createTempJSONFile({"radii": currentRadii, "slopes": currentSlopes}, "tempParams_")

    with multiprocessing.Pool() as pool:
        results = pool.starmap(processFunction, ((pdbid, currParamsFilename) for pdbid in pdbids))

    diffs = {atomType: [] for atomType in currentRadii}
    slopes = {atomType: [] for atomType in currentSlopes}
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

    diffs = { atomType:((analyser.medians.loc[atomType]['correctedDensity'] - analyser.chainMedian) / analyser.chainMedian) for atomType in radii if atomType in analyser.medians.index }
    newSlopes = { atomType:analyser.medians.loc[atomType]['slopes'] for atomType in slopes if atomType in analyser.medians.index }

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

