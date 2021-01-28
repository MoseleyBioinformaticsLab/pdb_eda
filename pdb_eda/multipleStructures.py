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
    --atom                          Aggregate and print results by atom
    --residue                       Aggregate and print results by residue
    --chain                         Aggregate and print results by chain
    --green                         Calculate and print results of all green blobs (positive difference electron density)
    --red                           Calculate and print results of all red blobs (negative difference electron density)
    --all                           Calculate and print results of both green and red blobs (positive and negative difference electron density)
    --stats                         If set true, return green or red blobs' statistics instead of blob object lists.
    --out-format=<format>           Output file format, available formats: csv, json [default: json].
    --symmetry-atoms                Calculate and print results of all symmetry atoms. (Only available in json format)
"""

import docopt
import os
import time
import sys

import json
import multiprocessing
import tempfile

from . import densityAnalysis
from . import __version__

## Final set of radii derived from optimization
defaultParamsPath = os.path.join(os.path.dirname(__file__), 'conf/optimized_radii_slope_param.json')
with open(defaultParamsPath, 'r') as fh:
    params = json.load(fh)

radiiGlobal = params['radii']
slopesGlobal = params['slopes']


def main():
    args = docopt.docopt(__doc__, version=__version__)
    if args["--help"]:
        print(__doc__)
        exit(0)

    pdbidFile = args['<pdbid-file>']
    resultFile = args['<out-file>']

    paramsPath = os.path.join(os.path.dirname(__file__), args["--radii-param"])
    if paramsPath != defaultParamsPath:
        with open(paramsPath, 'r') as fh:
            params = json.load(fh)
            radiiGlobal = params['radii']
            slopesGlobal = params['slopes']

    pdbids = []
    with open(pdbidFile, "r") as fileHandleIn:
        for pdbid in fileHandleIn:
            pdbids.append(pdbid[0:4])

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

    with open(resultFile, "w") as outfile:
        if args["--out-format"] == 'csv':
            statsHeaders = ['pdbid', 'chainMedian', 'voxelVolume', 'f000', 'chainNvoxel', 'chainTotalE', 'densityMean', 'diffDensityMean', 'resolution', 'spaceGroup']
            print(statsHeaders + sorted(radiiGlobal), sep=',', file=outfile)
            for result in fullResults.values():
                stats = [result["stats"][header] for header in statsHeader]
                diffs = [result["diffs"][atomType] for atomType in sorted(radiiGlobal)]
                print(stats + diffs, sep=',', file=outfile)
        else:
            json.dump(fullResults,outfile)


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

    analyser.aggregateCloud(radiiGlobal, slopesGlobal)
    analyser.estimateF000()
    if not analyser.chainMedian:
        return 0

    diffs = { atomType:((analyser.medians.loc[atomType]['correctedDensity'] - analyser.chainMedian) / analyser.chainMedian) if atomType in analyser.medians.index else 0 for atomType in sorted(radiiGlobal) }
    stats = { 'chainMedian' : analyser.chainMedian, 'voxelVolume' : analyser.densityObj.header.unitVolume, 'f000' : analyser.f000, 'chainNvoxel' : analyser.chainNvoxel, 'chainTotalE' : analyser.chainTotalE,
        'densityMean' : analyser.densityObj.header.densityMean, 'diffDensityMean' : analyser.diffDensityObj.header.densityMean, 'resolution' : analyser.pdbObj.header.resolution, 'spaceGroup' : analyser.pdbObj.header.spaceGroup }

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



