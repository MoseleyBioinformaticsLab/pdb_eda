# !/usr/bin/python3

import os
import time

import json
import multiprocessing
import tempfile

from . import densityAnalysis

## Final set of radii derived from optimization
radiiParamPath = os.path.join(os.path.dirname(__file__), 'conf/optimized_radii_slope_param.json')
with open(radiiParamPath, 'r') as fh:
    radiiParams = json.load(fh)

radiiDefault = radiiParams['radii']
slopesDefault = radiiParams['slopes']

def processFunction(pdbid, radii, slopes):
    analyser = densityAnalysis.fromPDBid(pdbid)

    if not analyser:
            return 0 

    analyser.aggregateCloud(radii, slopes)
    analyser.estimateF000()
    if not analyser.chainMedian:
        return 0

    diffs = []
    for atomType in sorted(radiiDefault):
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
    time.sleep(0.01) # sleep 0.01 seconds to prevent the same worker from calling this twice.
    return globalTempFile.name

def closeTemporaryFile(temp):
    filename = globalTempFile.name
    globalTempFile.close()
    # Sleep 0.01 seconds to prevent the same worker from calling this twice.
    # May need to increase this to 1 second or even 10 seconds to make this work.
    time.sleep(0.01)
    return filename


def main(args):
    pdbidFile = args['<pdbid-file>']
    resultFile = args['<out-file>']

    paramsPath = os.path.join(os.path.dirname(__file__), args["--radii-param"])
    with open(paramsPath, 'r') as fh:
        params = json.load(fh)
    radii = params['radii']
    slopes = params['slopes']

    pdbids = []
    with open(pdbidFile, "r") as fileHandleIn:
        for pdbid in fileHandleIn:
            pdbids.append(pdbid[0:4])

    with multiprocessing.Pool() as pool:
        result_filenames = pool.map(openTemporaryFile, (1 for i in range(multiprocessing.cpu_count())))
        results = pool.starmap(processFunction, ((pdbid, radii, slopes) for pdbid in pdbids))
        closed_filenames = pool.map(closeTemporaryFile, (1 for i in range(multiprocessing.cpu_count())))



    unclosed_filenames = set(result_filenames) - set(closed_filenames)
    if unclosed_filenames: # check whether all result files were closed.  Increase sleep time in closeTemporary File to prevent this.
        print("Unclosed Files: ", unclosed_filenames)

    with open(resultFile, "w") as outfile:
        print(*['pdbid', 'chainMedian', 'voxelVolume', 'f000', 'chainNvoxel', 'chainTotalE', 'densityMean', 'diffDensityMean', 'resolution', 'spaceGroup'] + sorted(radiiDefault), sep=',', file =outfile)
        for f in result_filenames:
            if f:
                with open(f, "r") as infile:
                    outfile.write(infile.read())
                os.remove(f) # remove the file


    #fileHandle = open(resultFile, 'w')
    #for result in results:
    #    if result:
    #        print(*result[1], *result[0], sep=", ", file=fileHandle)
