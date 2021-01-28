"""
pdb_eda multiple structure analysis mode command-line interface

Usage:
    pdb_eda multiple -h | --help
    pdb_eda multiple <pdbid-file> <out-file> [--radii-param=<paramfile>]

Options:
    -h, --help                      Show this screen.
    <pdbid>                         The PDB id
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
    --out-format=<format>           Onput file format, available formats: csv, json [default: json].
    --symmetry-atoms                Calculate and print results of all symmetry atoms. (Only available in jason format)
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
        result_filenames = pool.map(openTemporaryFile, (1 for i in range(multiprocessing.cpu_count())))
        results = pool.map(processFunction, pdbids)
        closed_filenames = pool.map(closeTemporaryFile, (1 for i in range(multiprocessing.cpu_count())))


    unclosed_filenames = set(result_filenames) - set(closed_filenames)
    if unclosed_filenames: # check whether all result files were closed.  Increase sleep time in closeTemporary File to prevent this.
        print("Error: unclosed intermediate files: ", unclosed_filenames, ". Results may not be complete.",file=sys.stderr)

    with open(resultFile, "w") as outfile:
        print(*['pdbid', 'chainMedian', 'voxelVolume', 'f000', 'chainNvoxel', 'chainTotalE', 'densityMean', 'diffDensityMean', 'resolution', 'spaceGroup'] + sorted(radiiDefault), sep=',', file=outfile)
        for filename in result_filenames:
            if filename:
                try:
                    with open(filename, "r") as infile:
                        outfile.write(infile.read())
                    os.remove(filename) # remove the file
                except:
                    print("Error: intermediate results file\"",filename,"\" was not parsable.  Results are not complete.",file=sys.stderr)


def processFunction(pdbid):
    analyser = densityAnalysis.fromPDBid(pdbid)

    if not analyser:
        return 0

    analyser.aggregateCloud(radiiGlobal, slopesGlobal)
    analyser.estimateF000()
    if not analyser.chainMedian:
        return 0

    diffs = []
    for atomType in sorted(radiiGlobal):
        diff = (analyser.medians.loc[atomType]['correctedDensity'] - analyser.chainMedian) / analyser.chainMedian if atomType in analyser.medians.index else 0
        diffs.append(diff)
    stats = [analyser.pdbid, analyser.chainMedian, analyser.densityObj.header.unitVolume, analyser.f000, analyser.chainNvoxel, analyser.chainTotalE,
             analyser.densityObj.header.densityMean, analyser.diffDensityObj.header.densityMean, analyser.pdbObj.header.resolution, analyser.pdbObj.header.spaceGroup]

    globalTempFile.write("%s\n" % ','.join([str(i) for i in stats + diffs]))
    return globalTempFile.name

def openTemporaryFile(temp):
    dirname = os.getcwd()
    global globalTempFile
    globalTempFile = tempfile.NamedTemporaryFile(mode='w', buffering=1, dir=dirname, prefix="tempPDB_", delete=False)
    time.sleep(0.1) # sleep 0.1 seconds to prevent the same worker from calling this twice.
    return globalTempFile.name

def closeTemporaryFile(temp):
    filename = globalTempFile.name
    globalTempFile.close()
    # Sleep 0.1 seconds to prevent the same worker from calling this twice.
    # May need to increase this to 1 second or even 10 seconds to make this work.
    time.sleep(0.1)
    return filename


