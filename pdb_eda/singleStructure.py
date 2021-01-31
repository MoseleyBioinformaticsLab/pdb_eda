# !/usr/bin/python3
"""
pdb_eda single structure analysis mode command-line interface
  Analyzes a single pdb entry.

Usage:
    pdb_eda single -h | --help
    pdb_eda single <pdbid> <out-file> map (--density | --diff-density)
    pdb_eda single <pdbid> <out-file> density (--atom | --residue | --chain) [--out-format=<format>] [--radii-param=<paramfile>]
    pdb_eda single <pdbid> <out-file> difference --blob [--green | --red] [--out-format=<format>] [--radii-param=<paramfile>]
    pdb_eda single <pdbid> <out-file> difference (--atom | --residue) [--green | --red] [--type=<type>] [--radius=<radius>] [--num-sd=<num-sd>] [--out-format=<format>] [--radii-param=<paramfile>]
    pdb_eda single <pdbid> <out-file> symmetry-atoms [--radii-param=<paramfile>]

Options:
    -h, --help                      Show this screen.
    <pdbid>                         The PDB ID to download and analyze.
    <out-file>                      Output filename. "-" will write to standard output.
    --radii-param=<paramfile>       Radii parameters. [default: conf/optimized_radii_slope_param.json]
    --atom                          Aggregate and print results by atom.
    --residue                       Aggregate and print results by residue.
    --chain                         Aggregate and print results by chain.
    --green                         Calculate and print results only for green blobs (positive difference electron density).
    --red                           Calculate and print results only for red blobs (negative difference electron density).
    --radius=<radius>               Radius (in angstroms) around atom or residue to calculate significant discrepancy. [default: 3.5]
    --num-sd=<num-sd>               Number of standard deviation units to use as a significant discrepancy cutoff. [default: 3]
    --type=<type>                   Residue type or atom type to filter by.
    --out-format=<format>           Output file format, available formats: csv, json [default: json].
    --symmetry-atoms                Calculate and print results of all symmetry atoms. (Only available in json format).
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
    print(args)
    exit(0)
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
    elif args["density"]:
        analyzer.aggregateCloud(radii, slopes, atomL=True, residueL=True, chainL=True)
        jsonType = "json"
        if args["--atom"]:
            headerList = list(map(str,list(analyzer.atomList) + ['density_electron_ratio']))
            result = [ list(item) + [analyzer.chainMedian] for item in analyzer.atomList.values.tolist()]
        elif args["--residue"]:
            headerList = ['chain', 'residue_number', 'residue_name', 'local_density_electron_ratio', 'num_voxels', 'electrons', 'volume', 'density_electron_ratio']
            result = [ list(item) + [analyzer.chainMedian] for item in analyzer.residueList]
        elif args["--chain"]:
            headerList = ['chain', 'residue_number', 'residue_name', 'local_density_electron_ratio', 'num_voxels', 'electrons', 'volumne', 'density_electron_ratio']
            result = [ list(item) + [analyzer.chainMedian] for item in analyzer.chainList]
    elif args["difference"]:
        jsonType = "json"
        if args["--atom"]:
            pass
        elif args["--residue"]:
            pass
        else:
            headerList = ['distance_to_atom', 'sign', 'electrons_of_discrepancy', 'num_voxels', 'volume', 'chain', 'residue_number', 'residue_name', 'atom_name', 'atom_symmetry', 'atom_xyz', 'centroid_xyz']
            result = analyzer.calcAtomBlobDists(radii, slopes)
            for blobInfo in result:
                blobInfo[10] = [ float(val) for val in blobInfo[10] ]
                blobInfo[11] = [ float(val) for val in blobInfo[11] ]
    elif args["symmetry-atoms"]:
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
