# !/usr/bin/python3
"""
pdb_eda single structure analysis mode command-line interface
  Analyzes a single pdb entry.

Usage:
    pdb_eda single -h | --help
    pdb_eda single <pdbid> <out-file> [--density-map | --diff-density-map]
    pdb_eda single <pdbid> <out-file> [--radii-param=<paramfile>] [--atom | --residue | --chain] [--out-format=<format>]
    pdb_eda single <pdbid> <out-file> [--radii-param=<paramfile>] [--green | --red | --all] [--stats] [--out-format=<format>]
    pdb_eda single <pdbid> <out-file> [--radii-param=<paramfile>] [--symmetry-atoms]

Options:
    -h, --help                      Show this screen.
    <pdbid>                         The PDB ID to download and analyze.
    <out-file>                      Output filename. "-" will write to standard output.
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
import sys
import json
import jsonpickle
from . import densityAnalysis
from . import __version__


def main():
    args = docopt.docopt(__doc__, version=__version__)
    if args["--help"]:
        print(__doc__)
        exit(0)

    paramsPath = os.path.join(os.path.dirname(__file__), args["--radii-param"])
    with open(paramsPath, 'r') as fh:
        params = json.load(fh)
    radii = params['radii']
    slopes = params['slopes']

    analyzer = densityAnalysis.fromPDBid(args["<pdbid>"])
    if not analyzer:
        sys.exit("Error: Unable to parse or download PDB entry or associated ccp4 file.")

    jsonType = "jsonpickle"
    if args["--density-map"]:
        result = analyzer.densityObj
    elif args["--diff-density-map"]:
        result = analyzer.diffDensityObj
    elif args["--atom"] or args["--residue"] or args["--chain"]:
        analyzer.aggregateCloud(radii, slopes, atomL=True, residueL=True, chainL=True)
        jsonType = "json"
        if args["--atom"]:
            headerList = list(map(str,list(analyzer.atomList) + ['chainMedianRatio']))
            result = [ list(item) + [analyzer.chainMedian] for item in analyzer.atomList.values.tolist()]
        elif args["--residue"]:
            headerList = ['chain', 'resNum', 'resName', 'density_electron_ratio', 'numVoxels', 'electrons', 'volume', 'chainMedianRatio']
            result = [ list(item) + [analyzer.chainMedian] for item in analyzer.residueList]
        elif args["--chain"]:
            headerList = ['chain', 'resNum', 'resName', 'density_electron_ratio', 'numVoxels', 'electrons', 'volumne', 'chainMedianRatio']
            result = [ list(item) + [analyzer.chainMedian] for item in analyzer.chainList]
    elif args["--green"] or args["--red"] or args["--all"]:
        if args["--stats"]:
            jsonType = "json"
            headerList = ['distance', 'sign', 'electrons_of_discrepancy', 'numVoxels', 'volume', 'chain', 'resNum', 'resName', 'atomName', 'atomSymmetry', 'atomXYZ', 'blobCentroid']
            result = analyzer.calcAtomBlobDists(radii, slopes)
            for blobInfo in result:
                blobInfo[10] = [ float(val) for val in blobInfo[10] ]
                blobInfo[11] = [ float(val) for val in blobInfo[11] ]
        else:
            analyzer.getBlobList()
            if args["--green"]:
                result = analyzer.greenBlobList
            if args["--red"]:
                result = analyzer.redBlobList
            if args["--all"]:
                result = analyzer.greenBlobList + analyzer.redBlobList
    elif args["--symmetry-atoms"]:
        analyzer.calcSymmetryAtoms()
        result = analyzer.symmetryAtoms

    with open(args["<out-file>"], 'w') if args["<out-file>"] != "-" else sys.stdout as outFile:
        if args["--out-format"] == 'csv':
            csvResults = [','.join(map(str, row)) for row in [headerList] + result]
            print(*csvResults, sep='\n', file=outFile)
        elif jsonType == "jsonpickle":
            jsonText = jsonpickle.encode(result)
            outFile.write(jsonText)
        else:
            jsonResults = [ dict(zip(headerList, row)) for row in result ]
            print(json.dumps(jsonResults, indent=2, sort_keys=True), file=outFile)
