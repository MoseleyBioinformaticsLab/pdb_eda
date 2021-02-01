# !/usr/bin/python3
"""
pdb_eda single structure analysis mode command-line interface
  Analyzes a single pdb entry.

Usage:
    pdb_eda single -h | --help
    pdb_eda single <pdbid> <out-file> map (--density | --diff-density)
    pdb_eda single <pdbid> <out-file> density (--atom | --residue | --chain) [--out-format=<format>] [--params=<params-file>] [--include-pdbid]
    pdb_eda single <pdbid> <out-file> difference --blob [--out-format=<format>] [--params=<params-file>] [--include-pdbid]
    pdb_eda single <pdbid> <out-file> difference (--atom | --residue) [--green | --red] [--type=<type>] [--radius=<radius>] [--num-sd=<num-sd>] [--out-format=<format>] [--params=<params-file>] [--include-pdbid]
    pdb_eda single <pdbid> <out-file> symmetry-atoms [--radii-param=<paramfile>]

Options:
    -h, --help                      Show this screen.
    <pdbid>                         The PDB ID to download and analyze.
    <out-file>                      Output filename. "-" will write to standard output.
    --params=<params-file>          Overriding parameters file that includes radii, slopes, etc. [default: ]
    --include-pdbid                 Include PDB ID at the beginning of each result.
    --atom                          Aggregate and print results by atom.
    --residue                       Aggregate and print results by residue.
    --chain                         Aggregate and print results by chain.
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

defaultParamsFilepath = os.path.join(os.path.dirname(__file__), 'conf/optimized_params.json')


def main():
    args = docopt.docopt(__doc__, version=__version__)
    args["--radius"] = float(args["--radius"])
    args["--num-sd"] = float(args["--num-sd"])

    if args["--help"]:
        print(__doc__)
        exit(0)

    paramsFilepath = args["--params"] if args["--params"] else defaultParamsFilepath
    try:
        with open(paramsFilepath, 'r') as paramsFile:
            params = json.load(paramsFile)
    except:
        sys.exit(str("Error: params file \"") + paramsFilepath + "\" does not exist or is not parsable.")

    analyzer = densityAnalysis.fromPDBid(args["<pdbid>"])
    if not analyzer:
        sys.exit("Error: Unable to parse or download PDB entry or associated ccp4 file.")

    jsonType = "jsonpickle"
    if args["--density"]:
        result = analyzer.densityObj
    elif args["--diff-density"]:
        result = analyzer.diffDensityObj
    elif args["density"]:
        analyzer.aggregateCloud(params, atomL=True, residueL=True, chainL=True)
        jsonType = "json"
        if args["--atom"]:
            headerList = list(map(str,list(analyzer.atomList) + ['density_electron_ratio']))
            result = [ list(item) + [analyzer.chainMedian] for item in analyzer.atomList.values.tolist()]
        elif args["--residue"]:
            headerList = densityAnalysis.DensityAnalysis.atomListHeader
            result = [ list(item) + [analyzer.chainMedian] for item in analyzer.residueList]
        elif args["--chain"]:
            headerList = densityAnalysis.DensityAnalysis.chainListHeader
            result = [ list(item) + [analyzer.chainMedian] for item in analyzer.chainList]
    elif args["difference"]:
        jsonType = "json"
        if args["--atom"]:
            headerList = densityAnalysis.DensityAnalysis.atomRegionDiscrepancyHeader
            result = analyzer.calculateAtomRegionDiscrepancies(args["--radius"], args["--num-sd"], args["--type"], params)
        elif args["--residue"]:
            headerList = densityAnalysis.DensityAnalysis.residueRegionDiscrepancyHeader
            result = analyzer.calculateResidueRegionDiscrepancies(args["--radius"], args["--num-sd"], args["--type"], params)
        else:
            headerList = densityAnalysis.DensityAnalysis.atomBlobDistanceHeader
            result = analyzer.calcAtomBlobDists(params)
            for blobInfo in result:
                blobInfo[10] = [ float(val) for val in blobInfo[10] ]
                blobInfo[11] = [ float(val) for val in blobInfo[11] ]
    elif args["symmetry-atoms"]:
        result = analyzer.symmetryAtoms

    if args["--include-pdbid"]:
        headerList = ["pdbid"] + headerList
        result = [ [analyzer.pdbid] + row for row in result ]

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
