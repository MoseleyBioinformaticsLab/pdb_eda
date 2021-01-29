"""
pdb_eda multiple structure analysis mode command-line interface
  Analyzes multiple pdb entries.

Usage:
    pdb_eda multiple -h | --help
    pdb_eda multiple <pdbid-file> <out-file> [--radii-param=<paramfile>] [--out-format=<format>]

Options:
    -h, --help                      Show this screen.
    <out-file>                      Output file name
    <pdbid-file>                    File name that contains the pdb ids
    --radii-param=<paramfile>       Radii parameters. [default: conf/optimized_radii_slope_param.json]
    --out-format=<format>           Output file format, available formats: csv, json [default: json].
"""

import docopt
import os
import time
import sys

import json
import csv
import multiprocessing
import tempfile

from . import densityAnalysis
from . import __version__

## Final set of radii derived from optimization
defaultParamsFilepath = os.path.join(os.path.dirname(__file__), 'conf/optimized_radii_slope_param.json')
with open(defaultParamsFilepath, 'r') as fh:
    params = json.load(fh)

radiiGlobal = params['radii']
slopesGlobal = params['slopes']
globalArgs = {}

def main():
    global globalArgs
    globalArgs = docopt.docopt(__doc__, version=__version__)
    if globalArgs["--help"]:
        print(__doc__)
        exit(0)

    paramsFilepath = os.path.join(os.path.dirname(__file__), globalArgs["--radii-param"])
    try:
        if paramsFilepath != defaultParamsFilepath:
            with open(paramsFilepath, 'r') as paramsFile:
                params = json.load(paramsFile)
                global radiiGlobal
                radiiGlobal = params['radii']
                global slopesGlobal
                slopesGlobal = params['slopes']
    except:
        sys.exit(str("Error: params file \"") + paramsFilepath + "\" does not exist or is not parsable.")

    try:
        pdbids = []
        with open(globalArgs["<pdbid-file>"], "r") as textFile:
            for pdbid in textFile:
                pdbids.append(pdbid[0:4])
    except:
        sys.exit(str("Error: PDB IDs file \"") + globalArgs["<pdbid-file>"] + "\" does not exist or is not parsable.")

    with multiprocessing.Pool() as pool:
        results = pool.map(processFunction, pdbids)

    fullResults = {}
    for resultFilename in results:
        if resultFilename:
            try:
                with open(resultFilename, 'r') as jsonFile:
                    result = json.load(jsonFile)
                    fullResults[result["pdbid"]] = result
                os.remove(resultFilename)
            except:
                pass

    if globalArgs["--out-format"] == 'csv':
        statsHeaders = ['chainMedian', 'voxelVolume', 'f000', 'chainNvoxel', 'chainTotalE', 'densityMean', 'diffDensityMean', 'resolution', 'spaceGroup', 'numAtomsAnalyzed', 'numResiduesAnalyzed', 'numChainsAnalyzed']
        with open(globalArgs['<out-file>'], "w", newline='') as csvFile:
            writer = csv.writer(csvFile)
            writer.writerow(['pdbid'] + statsHeaders + sorted(radiiGlobal))
            for result in fullResults.values():
                stats = [result["stats"][header] for header in statsHeaders]
                diffs = [result["diffs"][atomType] for atomType in sorted(radiiGlobal)]
                writer.writerow([result['pdbid']] + stats + diffs)
    else:
        with open(globalArgs['<out-file>'], "w") as jsonFile:
            json.dump(fullResults, jsonFile)


def processFunction(pdbid):
    """Process function to analyze a single pdb entry.

    :param :py:class:`str` pdbid: pdbid for entry to download and analyze.
    :return: resultFilename
    :rtype: :py:class:`str`
    """
    startTime = time.process_time()

    analyser = densityAnalysis.fromPDBid(pdbid)

    if not analyser:
        return 0

    analyser.aggregateCloud(radiiGlobal, slopesGlobal, atomL=True, residueL=True, chainL=True)
    analyser.estimateF000()
    if not analyser.chainMedian:
        return 0

    diffs = { atomType:((analyser.medians.loc[atomType]['correctedDensity'] - analyser.chainMedian) / analyser.chainMedian) if atomType in analyser.medians.index else 0 for atomType in sorted(radiiGlobal) }
    stats = { 'chainMedian' : analyser.chainMedian, 'voxelVolume' : analyser.densityObj.header.unitVolume, 'f000' : analyser.f000, 'chainNvoxel' : analyser.chainNvoxel, 'chainTotalE' : analyser.chainTotalE,
        'densityMean' : analyser.densityObj.header.densityMean, 'diffDensityMean' : analyser.diffDensityObj.header.densityMean, 'resolution' : analyser.pdbObj.header.resolution,
        'spaceGroup' : analyser.pdbObj.header.spaceGroup, 'numAtomsAnalyzed' : len(analyser.atomList.index), 'numResiduesAnalyzed' : len(analyser.residueList), 'numChainsAnalyzed' : len(analyser.chainList)  }

    elapsedTime = time.process_time() - startTime
    resultFilename = createTempJSONFile({ "pdbid" : analyser.pdbid, "diffs" : diffs, "stats" : stats, "execution_time" : elapsedTime }, "tempResults_")
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



