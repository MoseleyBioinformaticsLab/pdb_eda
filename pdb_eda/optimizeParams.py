#!/usr/bin/python3
"""
The radii and slope parameter optimization mode command-line interface
  Optimizes radii and b-factor slopes using a given set of PDB IDs.
  A simple steepest decent optimization approach is utilized.
  This approach is justified by testing and the use of median differences that smooths the error surface.

Usage:
    pdb_eda optimize -h | --help
    pdb_eda optimize <start-params-file> <pdbid-file> <log-file> <out-params-file> [options]
    pdb_eda optimize <params-file1> <params-file2> --compare
    pdb_eda optimize <start-params-file> <out-params-file> --finalize

Options:
    -h, --help                          Show this screen.
    --ignore                            Ignore the "optimize" atom type limit in the parameter file.
    --reverse                           Reverse the "optimize" atom type limit.
    --sample=<sample-size>              Use a random sample of PDB ids for optimization. [default: 0]
    --max=<max-radius-change>           Maximum to change the radius at each incremental optimization. [default: 0.2]
    --min=<min-radius-change>           Minimum to change the radius at each incremental optimization. [default: 0.001]
    --radius=<start-radius>             Starting radius for the starting atom-type. [default: 0]
    --start=<start-atom-type>           Starting atom type. [default: ]
    --stop=<fractional-difference>      Max penalty fraction allowed for stopping the optimization. [default: 0]
    --unweighted                        Pick atom type to optimize next without weighting based on occurrence across PDB entries.
    --penalty-weight=<inverse-weight>   Inverse penalty weight for the atom-type specific overlap completeness. [default: 3.0]
    --compare                           Compare two parameter files.
    --finalize                          Finalize the parameter file for general use.
    --testing                           Run only a single process for testing purposes.


The optimization uses the absolute value of a penalty function to determine which atom-type radius to optimize next.
The penalty function is: medianDensityElectronRatioDifferences[atomType] + (overlapCompleteness[atomType] - maxOverlapCompleteness)/inversePenaltyWeight
    overlapCompleteness - fraction of atoms of a specific atom-type where the density cloud around the atom overlaps with the density clouds of bonded atoms.
    maxOverlapCompleteness - maximum of all atom-type-specific overlap completeness values.
    The median density-electron ratio differences tends to optimize smaller radii, while overlap completeness tends to optimize larger radii.
    The two terms of the penalty function balance two opposing criteria that need to be maximized to create the best estimate of the density-electron ratio.
    Also negative penalty function values indicate that a positive change in radius is needed, while a positive penalty function values indicate that a negative change in radius is needed. 

This mode is often run multiple times using the output parameter file generated in one cycle as the starting parameter file in the next cycle.
Typically, it is good to start with the following series of cycles.  You can start with larger sample sizes if you have more than 20 CPU cores available.
The total optimization time is typically between 5000 to 10000 CPU hours.
1) --sample=100 --max=0.1 --min=0.001
2) --sample=100 --max=0.1 --min=0.001 --unweighted
5) (use ~1000 pdbids) --max=0.1 --min=0.001
6) (use ~1000 pdbids) --max=0.1 --min=0.001 --unweighted
If you need to kill the run, don't worry, there is a temp output parameter file with the last improvement.
Also, it is good to add prior amino acid optimized atom types to additional atom types being optimized.  See generate mode --params=<params-file> option for details.
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

    if args["--compare"]:
        try:
            with open(args["<params-file1>"], 'r') as jsonFile:
                params1 = json.load(jsonFile)
        except:
            RuntimeError(str("Error: params file \"") + args["<params-file1>"] + "\" does not exist or is not parsable.")

        try:
            with open(args["<params-file2>"], 'r') as jsonFile:
                params2 = json.load(jsonFile)
        except:
            RuntimeError(str("Error: params file \"") + args["<params-file2>"] + "\" does not exist or is not parsable.")


        atomTypes = set(params1["radii"].keys()).union(params2["radii"].keys())
        radiusDifferences = {atomType:params1["radii"][atomType] - params2["radii"][atomType] for atomType in atomTypes
                             if atomType in params1["radii"] and not np.isnan(params1["radii"][atomType]) and atomType in params2["radii"] and not np.isnan(params2["radii"][atomType])}
        maxRadiusDiffAtomType = max(radiusDifferences, key=lambda y: abs(radiusDifferences[y]))
        meanRadiusDifference = np.nanmean(list(radiusDifferences.values()))
        StDRadiusDifferences = np.nanstd(list(radiusDifferences.values()))
        print("Radii Comparison:",args["<params-file1>"],"vs", args["<params-file2>"])
        print("Max Radius Difference:",radiusDifferences[maxRadiusDiffAtomType],"for",maxRadiusDiffAtomType,", leaving_atom =", maxRadiusDiffAtomType in params1["leaving_atoms"])
        print("Mean (Std) Radius Differences:", meanRadiusDifference,str("(")+str(StDRadiusDifferences)+")")

        NanAtomTypes = [atomType for (atomType,radius) in params1["radii"].items() if np.isnan(radius)]
        if NanAtomTypes:
            print("AtomTypes in",args["<params-file1>"], "with NaN radius:", ", ".join(NanAtomTypes))
        NanAtomTypes = [atomType for (atomType, radius) in params2["radii"].items() if np.isnan(radius)]
        if NanAtomTypes:
            print("AtomTypes in", args["<params-file2>"], "with NaN radius:", ", ".join(NanAtomTypes))
        NanAtomTypes = [atomType for (atomType,slope) in params1["slopes"].items() if np.isnan(slope)]
        if NanAtomTypes:
            print("AtomTypes in",args["<params-file1>"], "with NaN slope:", ", ".join(NanAtomTypes))
        NanAtomTypes = [atomType for (atomType, slope) in params2["slopes"].items() if np.isnan(slope)]
        if NanAtomTypes:
            print("AtomTypes in", args["<params-file2>"], "with NaN slope:", ", ".join(NanAtomTypes))
    elif args["--finalize"]:
        try:
            with open(args["<start-params-file>"], 'r') as jsonFile:
                params = json.load(jsonFile)
        except:
            RuntimeError(str("Error: params file \"") + args["<start-params-file>"] + "\" does not exist or is not parsable.")

        if "optimize" in params:
            del params["optimize"]

        try:
            with open(args["<out-params-file>"], 'w') as jsonFile:
                print(json.dumps(params, indent=2, sort_keys=True), file=jsonFile)
        except:
            RuntimeError(str("Error: unable to create params file \"") + args["<out-params-file>"] + "\".")
    else:
        maxRadiusIncrement = float(args["--max"])
        radiusIncrement = maxRadiusIncrement
        minRadiusIncrement = float(args["--min"])
        stoppingFractionalDifference = float(args["--stop"])
        startingRadius = float(args["--radius"])
        sampleSize = int(args["--sample"])
        inversePenaltyWeight = float(args["--penalty-weight"])
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
            RuntimeError(str("Error: params file \"") + args["<start-params-file>"] + "\" does not exist or is not parsable.")

        if args["--reverse"] and atomTypes2Optimize:
            atomTypes2Optimize = { atomType for atomType in currentRadii.keys() if atomType not in atomTypes2Optimize }

        if args["--start"] != "" and args["--start"] not in currentRadii:
            RuntimeError(str("Error: starting atom \"") + args["--start"] + "\" is not valid.")

        try:
            pdbids = []
            with open(args["<pdbid-file>"], "r") as textFile:
                for pdbid in textFile:
                    pdbids.append(pdbid[0:4])
        except:
            RuntimeError(str("Error: PDB IDs file \"") + args["<pdbid-file>"] + "\" does not exist or is not parsable.")

        if sampleSize > 0:
            pdbids = random.sample(pdbids,sampleSize)

        with open(args["<log-file>"], 'w') as logFile:
            print(args, file=logFile)

            print("PDB IDs:",",".join(pdbids), file=logFile)
            print("Calculating start median differences: start-time=", str(datetime.datetime.now()))
            print("Calculating start median differences: start-time=", str(datetime.datetime.now()), file=logFile)

            (bestMedianDiffs, meanDiffs, overallStdDevDiffs, currentSlopes, sizes, overlapCompleteness) = calculateMedianDiffsSlopes(pdbids, params, args["--testing"], args["<pdbid-file>"]+".execution_times")
            currentSlopes = {**currentSlopes, **(params["slopes"])}
            maxOverlapCompleteness = max(overlapCompleteness.values())
            bestPenalties = {atomType:(bestMedianDiffs[atomType] + (overlapCompleteness[atomType] - maxOverlapCompleteness)/inversePenaltyWeight) for atomType in bestMedianDiffs}

            maxSize = max([sizes[atomType] for atomType in bestMedianDiffs.keys() if not atomTypes2Optimize or atomType in atomTypes2Optimize])
            print("Starting Radii Min-Max: [",min(currentRadii.values()),",",max(currentRadii.values()),"]",file=logFile)
            print("Max Absolute Weighted Median Diff:", max([abs(bestMedianDiffs[atomType] * sizes[atomType] / maxSize) for atomType in bestMedianDiffs.keys() if not atomTypes2Optimize or atomType in atomTypes2Optimize]),
                  ", Weighted Diff StdDev:", overallStdDevDiffs,
                  ", Max Size:", maxSize)
            print("Max Absolute Weighted Median Diff:", max([abs(bestMedianDiffs[atomType] * sizes[atomType] / maxSize) for atomType in bestMedianDiffs.keys() if not atomTypes2Optimize or atomType in atomTypes2Optimize]),
                  ", Weighted Diff StdDev:", overallStdDevDiffs,
                  ", Max Size:", maxSize, file=logFile)
            print("Max Absolute Median Diff:", max([abs(bestMedianDiffs[atomType]) for atomType in bestMedianDiffs.keys() if not atomTypes2Optimize or atomType in atomTypes2Optimize]),
                  ", Max Abs Diff Mean-Median:", max([abs(meanDiffs[atomType] - bestMedianDiffs[atomType]) for atomType in bestMedianDiffs.keys() if not atomTypes2Optimize or atomType in atomTypes2Optimize]),
                  ", Mean Abs Diff Mean-Median:", np.mean([abs(meanDiffs[atomType] - bestMedianDiffs[atomType]) for atomType in bestMedianDiffs.keys() if not atomTypes2Optimize or atomType in atomTypes2Optimize]))
            print("Max Absolute Median Diff:", max([abs(bestMedianDiffs[atomType]) for atomType in bestMedianDiffs.keys() if not atomTypes2Optimize or atomType in atomTypes2Optimize]),
                  ", Max Abs Diff Mean-Median:", max([abs(meanDiffs[atomType] - bestMedianDiffs[atomType]) for atomType in bestMedianDiffs.keys() if not atomTypes2Optimize or atomType in atomTypes2Optimize]),
                  ", Mean Abs Diff Mean-Median:", np.mean([abs(meanDiffs[atomType] - bestMedianDiffs[atomType]) for atomType in bestMedianDiffs.keys() if not atomTypes2Optimize or atomType in atomTypes2Optimize]), file=logFile)
            print("Max Absolute Weighted Penalty:", max([abs(bestPenalties[atomType] * sizes[atomType] / maxSize) for atomType in bestPenalties.keys() if not atomTypes2Optimize or atomType in atomTypes2Optimize]),
                  ", max overlap completeness=",maxOverlapCompleteness)
            print("Max Absolute Weighted Penalty:", max([abs(bestPenalties[atomType] * sizes[atomType] / maxSize) for atomType in bestPenalties.keys() if not atomTypes2Optimize or atomType in atomTypes2Optimize]),
                  ", max overlap completeness=",maxOverlapCompleteness, file=logFile)
            print("Overlap Completeness Min-Max: [", min(overlapCompleteness.values()),",",max(overlapCompleteness.values()),"]")
            print("Overlap Completeness Min-Max: [", min(overlapCompleteness.values()),",",max(overlapCompleteness.values()),"]", file=logFile)
            print("Radii:", currentRadii, file=logFile)
            print("Median Diffs:", bestMedianDiffs, file=logFile)
            print("Overlap Completeness:", overlapCompleteness, file=logFile)
            print("Penalties:", bestPenalties, file=logFile)


            testBestPenalties = {atomType: penalty for (atomType, penalty) in bestPenalties.items() if atomType in atomTypes2Optimize} if atomTypes2Optimize else bestPenalties
            if args["--unweighted"]:
                currentAtomType = max(testBestPenalties, key=lambda y: abs(testBestPenalties[y])) if not args["--start"]  else args["--start"]
            else:
                currentAtomType = max(testBestPenalties, key=lambda y: abs(testBestPenalties[y] * sizes[y])) if not args["--start"] else args["--start"]
            previousRadius = currentRadii[currentAtomType]

            if startingRadius > 0:
                previousDirection = currentRadii[currentAtomType] < startingRadius
                currentRadii[currentAtomType] = startingRadius
            else:
                currentRadii[currentAtomType] = currentRadii[currentAtomType] + radiusIncrement if bestPenalties[currentAtomType] < 0 else currentRadii[currentAtomType] - radiusIncrement
                previousDirection = bestPenalties[currentAtomType] < 0

            numAccepted = 0
            numRejected = 0
            estimatedRadiusIncrement = {atomType:0 for atomType in currentRadii.keys()}
            while True:
                print("Testing ", currentAtomType, ": starting radius=", previousRadius, ", new radius=", currentRadii[currentAtomType], ", current weighted penalty=", bestPenalties[currentAtomType] * sizes[currentAtomType] / maxSize,
                      ", current weighted median difference=",bestMedianDiffs[currentAtomType] * sizes[currentAtomType] / maxSize, str("(")+str(bestMedianDiffs[currentAtomType])+str(")"), ", size=",sizes[currentAtomType])
                print("Testing ", currentAtomType, ": starting radius=", previousRadius, ", new radius=", currentRadii[currentAtomType], ", current weighted penalty=", bestPenalties[currentAtomType] * sizes[currentAtomType] / maxSize,
                      ", current median difference=",bestMedianDiffs[currentAtomType], str("(")+str(bestMedianDiffs[currentAtomType])+str(")"), ", size=",sizes[currentAtomType], file=logFile)
                print("Calculating next  median differences: start-time=", str(datetime.datetime.now()), ", current increment=",radiusIncrement)
                print("Calculating next  median differences: start-time=", str(datetime.datetime.now()), ", current increment=",radiusIncrement, file=logFile)

                (medianDiffs, meanDiffs, overallStdDevDiffs, slopes, sizes, overlapCompleteness) = calculateMedianDiffsSlopes(pdbids, {**params, "radii" : currentRadii, "slopes" : currentSlopes }, args["--testing"],
                                                                                                         args["<pdbid-file>"]+".execution_times")
                maxOverlapCompleteness = max(overlapCompleteness.values())
                penalties = {atomType:(medianDiffs[atomType] + (overlapCompleteness[atomType] - maxOverlapCompleteness)/inversePenaltyWeight) for atomType in medianDiffs}

                maxSize = max([sizes[atomType] for atomType in medianDiffs.keys() if not atomTypes2Optimize or atomType in atomTypes2Optimize])
                print("Radii:", currentRadii, file=logFile)
                print("Median Diffs:", medianDiffs, file=logFile)
                print("Overlap Completeness:", overlapCompleteness, file=logFile)
                print("Penalties:", penalties, file=logFile)
                print("Slopes:", slopes, file=logFile)
                print("Max Absolute Weighted Median Diff:", max([abs(medianDiffs[atomType] * sizes[atomType] / maxSize) for atomType in medianDiffs.keys() if not atomTypes2Optimize or atomType in atomTypes2Optimize]),
                      ", Weighted Diff StdDev:", overallStdDevDiffs,
                      ", Max Size:", maxSize)
                print("Max Absolute Weighted Median Diff:", max([abs(medianDiffs[atomType] * sizes[atomType] / maxSize) for atomType in medianDiffs.keys() if not atomTypes2Optimize or atomType in atomTypes2Optimize]),
                      ", Weighted Diff StdDev:", overallStdDevDiffs,
                      ", Max Size:", maxSize, file=logFile)
                print("Max Absolute Median Diff:", max([abs(medianDiffs[atomType]) for atomType in medianDiffs.keys() if not atomTypes2Optimize or atomType in atomTypes2Optimize]),
                      ", Max Abs Diff Mean-Median:", max([abs(meanDiffs[atomType] - medianDiffs[atomType]) for atomType in medianDiffs.keys() if not atomTypes2Optimize or atomType in atomTypes2Optimize]),
                      ", Mean Abs Diff Mean-Median:", np.mean([abs(meanDiffs[atomType] - medianDiffs[atomType]) for atomType in medianDiffs.keys() if not atomTypes2Optimize or atomType in atomTypes2Optimize]))
                print("Max Absolute Median Diff:", max([abs(medianDiffs[atomType]) for atomType in medianDiffs.keys() if not atomTypes2Optimize or atomType in atomTypes2Optimize]),
                      ", Max Abs Diff Mean-Median:", max([abs(meanDiffs[atomType] - medianDiffs[atomType]) for atomType in medianDiffs.keys() if not atomTypes2Optimize or atomType in atomTypes2Optimize]),
                      ", Mean Abs Diff Mean-Median:", np.mean([abs(meanDiffs[atomType] - medianDiffs[atomType]) for atomType in medianDiffs.keys() if not atomTypes2Optimize or atomType in atomTypes2Optimize]), file=logFile)
                print("Max Absolute Weighted Penalty:", max([abs(penalties[atomType] * sizes[atomType] / maxSize) for atomType in penalties.keys() if not atomTypes2Optimize or atomType in atomTypes2Optimize]),
                      ", max overlap completeness=",maxOverlapCompleteness)
                print("Max Absolute Weighted Penalty:", max([abs(penalties[atomType] * sizes[atomType] / maxSize) for atomType in penalties.keys() if not atomTypes2Optimize or atomType in atomTypes2Optimize]),
                      ", max overlap completeness=",maxOverlapCompleteness, file=logFile)

                improved = False
                directionChangeByIncrement = (previousDirection != (penalties[currentAtomType] < 0)) and estimatedRadiusIncrement[currentAtomType] == 0
                if abs(penalties[currentAtomType]) <= abs(bestPenalties[currentAtomType]):
                    numAccepted += 1
                    if abs(penalties[currentAtomType]) < abs(bestPenalties[currentAtomType]):
                        estimatedRadiusIncrement[currentAtomType] = 0.9 * (currentRadii[currentAtomType] - previousRadius) * penalties[currentAtomType] / (bestPenalties[currentAtomType] - penalties[currentAtomType])
                        if abs(estimatedRadiusIncrement[currentAtomType]) < minRadiusIncrement:
                            estimatedRadiusIncrement[currentAtomType] = 0
                    else:
                        estimatedRadiusIncrement[currentAtomType] = 0
                    bestMedianDiffs = medianDiffs
                    bestPenalties = penalties
                    currentSlopes = {**slopes, **currentSlopes}
                    improved = True if abs(penalties[currentAtomType]) < abs(bestPenalties[currentAtomType]) else 2

                    print("Accepted", currentAtomType, ": new radius=", currentRadii[currentAtomType], ", current weighted penalty=", bestPenalties[currentAtomType] * sizes[currentAtomType] / maxSize,
                          ", current weighted median difference=",bestMedianDiffs[currentAtomType] * sizes[currentAtomType] / maxSize,
                          str("(")+str(bestMedianDiffs[currentAtomType])+str(")"), ", size=",sizes[currentAtomType])
                    print("Accepted", currentAtomType, ": new radius=", currentRadii[currentAtomType], ", current weighted penalty=", bestPenalties[currentAtomType] * sizes[currentAtomType] / maxSize,
                          ", current weighted median difference=",bestMedianDiffs[currentAtomType] * sizes[currentAtomType] / maxSize,
                          str("(")+str(bestMedianDiffs[currentAtomType])+str(")"), ", size=",sizes[currentAtomType], file=logFile)

                    try:
                        with open(args["<out-params-file>"] + ".temp", 'w') as jsonFile:
                            print(json.dumps({**params, "radii": currentRadii, "slopes": currentSlopes}, indent=2, sort_keys=True), file=jsonFile)
                    except:
                        RuntimeError(str("Error: unable to create temporary params file \"") + args["<out-params-file>"] + ".temp" + "\".")
                else:
                    numRejected += 1
                    estimatedRadiusIncrement[currentAtomType] = 0
                    print("Rejected", currentAtomType, ": new radius=", currentRadii[currentAtomType])
                    print("Rejected", currentAtomType, ": new radius=", currentRadii[currentAtomType], file=logFile)
                    currentRadii[currentAtomType] = previousRadius

                testBestPenalties = { atomType:diff for (atomType,diff) in bestPenalties.items() if atomType in atomTypes2Optimize } if atomTypes2Optimize else bestPenalties
                if args["--unweighted"]:
                    maxAtomType = max(testBestPenalties, key=lambda y: abs(testBestPenalties[y]))
                else:
                    maxAtomType = max(testBestPenalties, key=lambda y: abs(testBestPenalties[y]) * sizes[y])

                if stoppingFractionalDifference > 0 and max([abs(value * sizes[atomType] / maxSize) for atomType,value in testBestPenalties.items()]) < stoppingFractionalDifference:
                    break

                if maxAtomType == currentAtomType:
                    if not improved or previousDirection != (bestPenalties[currentAtomType] < 0):
                        if radiusIncrement == minRadiusIncrement:
                            break

                        radiusIncrement = radiusIncrement / 2.0
                        if radiusIncrement < minRadiusIncrement:
                            radiusIncrement = minRadiusIncrement
                    elif improved == 2:
                        radiusIncrement = radiusIncrement * 1.5
                        if radiusIncrement > maxRadiusIncrement:
                            radiusIncrement = maxRadiusIncrement

                elif directionChangeByIncrement:
                    radiusIncrement = radiusIncrement * 0.9
                    if radiusIncrement < minRadiusIncrement:
                        break

                currentAtomType = maxAtomType
                previousRadius = currentRadii[currentAtomType]
                if abs(estimatedRadiusIncrement[currentAtomType]) > 0:
                    currentRadii[currentAtomType] = currentRadii[currentAtomType] + estimatedRadiusIncrement[currentAtomType]
                else:
                    currentRadii[currentAtomType] = currentRadii[currentAtomType] + radiusIncrement if bestPenalties[currentAtomType] < 0 else currentRadii[currentAtomType] - radiusIncrement
                previousDirection = bestPenalties[currentAtomType] < 0
                gc.collect() # force garbage collection to decrease memory use.

            print("Final Radii:", currentRadii)
            print("Final Radii:", currentRadii, file=logFile)
            print("Final Radii Min-Max: [",min(currentRadii.values()),",",max(currentRadii.values()),"]")
            print("Final Radii Min-Max: [",min(currentRadii.values()),",",max(currentRadii.values()),"]",file=logFile)
            print("Num Accepted Changes=", numAccepted, ", Num Rejected Changes=", numRejected)
            print("Num Accepted Changes=", numAccepted, ", Num Rejected Changes=", numRejected, file=logFile)
            print("Max Absolute Weighted Median Diff:", max([abs(bestMedianDiffs[atomType] * sizes[atomType] / maxSize) for atomType in bestMedianDiffs.keys() if not atomTypes2Optimize or atomType in atomTypes2Optimize]))
            print("Max Absolute Weighted Median Diff:", max([abs(bestMedianDiffs[atomType] * sizes[atomType] / maxSize) for atomType in bestMedianDiffs.keys() if not atomTypes2Optimize or atomType in atomTypes2Optimize]), file=logFile)
            print("Max Absolute Weighted Penalty:", max([abs(testBestPenalties[atomType] * sizes[atomType] / maxSize) for atomType in testBestPenalties.keys()]))
            print("Max Absolute Weighted Penalty:", max([abs(testBestPenalties[atomType] * sizes[atomType] / maxSize) for atomType in testBestPenalties.keys()]), file=logFile)
            print("Overlap Completeness Min-Max: [", min(overlapCompleteness.values()),",",max(overlapCompleteness.values()),"]")
            print("Overlap Completeness Min-Max: [", min(overlapCompleteness.values()),",",max(overlapCompleteness.values()),"]", file=logFile)
            print("Optimization end-time=", str(datetime.datetime.now()))
            print("Optimization end-time=", str(datetime.datetime.now()), file=logFile)

            outParams = { **params, "radii" : currentRadii, "slopes" : currentSlopes }

        try:
            with open(args["<out-params-file>"], 'w') as jsonFile:
                print(json.dumps(outParams, indent=2, sort_keys=True), file=jsonFile)
        except:
            RuntimeError(str("Error: unable to create params file \"") + args["<out-params-file>"] + "\".")

def calculateMedianDiffsSlopes(pdbids, currentParams, testing=False, executionTimesFilename=None):
    """Calculates the median diffs and slopes across a list of pdb entries.

    :param pdbids: list of pdbids to process.
    :type pdbids: :py:class:`list`
    :param currentParams:  parameters.
    :type currentParams: :py:class:`dict`

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
    atomTypeOverlapCompleteness = {atomType: 0 for atomType in currentParams["radii"]}
    atomTypeOverlapIncompleteness = {atomType: 0 for atomType in currentParams["radii"]}

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
                    for atomType, count in result['atomtype_overlap_completeness'].items():
                        atomTypeOverlapCompleteness[atomType] += count
                    for atomType, count in result['atomtype_overlap_incompleteness'].items():
                        atomTypeOverlapIncompleteness[atomType] += count
                os.remove(resultFilename)
            except:
                pass

    for atomType in atomTypeOverlapCompleteness.keys():
        if atomTypeOverlapCompleteness[atomType] > 0 or atomTypeOverlapIncompleteness[atomType] > 0:
            atomTypeOverlapCompleteness[atomType] = atomTypeOverlapCompleteness[atomType] / (atomTypeOverlapCompleteness[atomType] + atomTypeOverlapIncompleteness[atomType])
        else:
            atomTypeOverlapCompleteness[atomType] = 1 # This makes the overlap penalty zero.

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
    medianSlopes = {key:value for (key,value) in medianSlopes.items() if not np.isnan(value)} # remove nan values.

    return (medianDiffs, meanDiffs, overallStdDevDiffs, medianSlopes, sizeDiffs, atomTypeOverlapCompleteness)

def processFunction(pdbid, paramsFilepath):
    """Process function to analyze a single pdb entry.

    :param pdbid: pdbid for entry to download and analyze.
    :type pdbid: :py:class:`str`
    :param paramsFilepath: filepath to the radii and slopes parameters.
    :type paramsFilepath: :py:class:`str`

    :return: resultFilename or 0
    :rtype: :py:class:`str`, :py:class:`int`
    """
    try:
        with open(paramsFilepath, 'r') as jsonFile:
            params = json.load(jsonFile)
        densityAnalysis.setGlobals(params)
    except:
        return 0

    startTime = time.process_time()

    analyzer = densityAnalysis.fromPDBid(pdbid)
    if not analyzer or not analyzer.densityElectronRatio:
        return 0

    diffs = {atomType:((analyzer.medians['corrected_density_electron_ratio'][atomType] - analyzer.densityElectronRatio) / analyzer.densityElectronRatio) for atomType in params["radii"]
             if atomType in analyzer.medians['corrected_density_electron_ratio'] and not np.isnan(analyzer.medians['corrected_density_electron_ratio'][atomType])}
    newSlopes = {atomType:analyzer.medians['slopes'][atomType] for atomType in params["slopes"] if atomType in analyzer.medians['slopes'] and not np.isnan(analyzer.medians['slopes'][atomType])}

    elapsedTime = time.process_time() - startTime
    resultFilename = fileUtils.createTempJSONFile({ "pdbid" : pdbid, "diffs" : diffs, "slopes" : newSlopes, 'resolution' : analyzer.pdbObj.header.resolution, "execution_time" : elapsedTime,
                                                    "atomtype_overlap_completeness" : analyzer.atomTypeOverlapCompleteness, "atomtype_overlap_incompleteness" : analyzer.atomTypeOverlapIncompleteness }, "tempResults_")

    # force garbage collection to decrease memory use.
    analyzer=0
    diffs=0
    newSlopes=0
    gc.collect()

    return resultFilename


