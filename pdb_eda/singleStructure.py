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
    <pdbid>                         The PDB id
    <out-file>                      Output file name
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
    pdbid = args["<pdbid>"]
    filename = args["<out-file>"]

    paramsPath = os.path.join(os.path.dirname(__file__), args["--radii-param"])
    with open(paramsPath, 'r') as fh:
        params = json.load(fh)
    radii = params['radii']
    slopes = params['slopes']

    analyser = densityAnalysis.fromPDBid(pdbid)
    if not analyser:
        sys.exit("Error: Unable to parse or download PDB entry or associated ccp4 file.")

    result = []
    if args["--density-map"]:
        result = analyser.densityObj
    elif args["--diff-density-map"]:
        result = analyser.diffDensityObj
    elif args["--atom"] or args["--residue"] or args["--chain"]:
        analyser.aggregateCloud(radii, slopes, atomL=True, residueL=True, chainL=True)
        if args["--atom"]:
            result.append(','.join(map(str, list(analyser.atomList) + ['chainMedianRatio'])))
            for item in analyser.atomList.values.tolist():
                result.append(','.join(map(str, item + [analyser.chainMedian])))
        elif args["--residue"]:
            result.append(','.join(['chain', 'resNum', 'resName', 'density_electron_ratio', 'volume', 'electrons', 'chainMedianRatio']))
            for item in analyser.residueList:
                result.append(','.join(map(str, item + [analyser.chainMedian])))
        elif args["--chain"]:
            result.append(','.join(['chain', 'resNum', 'resName', 'density_electron_ratio', 'volume', 'electrons', 'chainMedianRatio']))
            for item in analyser.chainList:
                result.append(','.join(map(str, item + [analyser.chainMedian])))
    elif args["--green"] or args["--red"] or args["--all"]:
        if args["--stats"]:
            for item in analyser.calcAtomBlobDists(radii, slopes):
                result.append(','.join(map(str, item)))
        else:
            analyser.getBlobList()
            if args["--green"]:
                result = analyser.greenBlobList
            if args["--red"]:
                result = analyser.redBlobList
            if args["--all"]:
                result = analyser.greenBlobList + analyser.redBlobList
    elif args["--symmetry-atoms"]:
        analyser.calcSymmetryAtoms()
        result = analyser.symmetryAtoms


    with open(filename, 'w') as outFile:
        if args["--out-format"] == 'csv':
            print(*result, sep='\n', file=outFile)
        else:
            result = jsonpickle.encode(result)
            outFile.write(result)
