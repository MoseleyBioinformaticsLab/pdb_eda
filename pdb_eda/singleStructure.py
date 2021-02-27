# !/usr/bin/python3
"""
pdb_eda single structure analysis mode command-line interface
  Analyzes a single pdb entry.

Usage:
    pdb_eda single -h | --help
    pdb_eda single <pdbid> <out-file> map (--density | --diff-density)
    pdb_eda single <pdbid> <out-file> density (--atom | --residue | --domain) [--out-format=<format>] [--params=<params-file>] [--include-pdbid]
    pdb_eda single <pdbid> <out-file> difference (--atom | --residue | --symmetry-atom) [--type=<type>] [--radius=<radius>] [--num-sd=<num-sd>] [--out-format=<format>] [--params=<params-file>] [--include-pdbid]
    pdb_eda single <pdbid> <out-file> blob [--green] [--red] [--num-sd=<num-sd>] [--out-format=<format>] [--params=<params-file>] [--include-pdbid]
    pdb_eda single <pdbid> <out-file> blob --blue [--num-sd=<num-sd>] [--out-format=<format>] [--params=<params-file>] [--include-pdbid]
    pdb_eda single <pdbid> <out-file> statistics [--out-format=<format>] (--atom | --residue) [--include-pdbid] [--print-validation]

Options:
    -h, --help                      Show this screen.
    <pdbid>                         The PDB ID to download and analyze.
    <out-file>                      Output filename. "-" will write to standard output.
    --params=<params-file>          Overriding parameters file that includes radii, slopes, etc. [default: ]
    --include-pdbid                 Include PDB ID at the beginning of each result.
    --density                       Output the density Fo-Fc map in jsonpickle format.
    --diff-density                  Output the difference density 2Fo-Fc map in jsonpickle format.
    --atom                          Calculate results for each atom.
    --residue                       Calculate results for each residue.
    --symmetry-atom                 Calculate results for each symmetry atom.
    --domain                        Calculate results for each domain.
    --green                         Calculate green (positive) difference map blobs.
    --red                           Calculate red (negative) difference map blobs.
    --blue                          Calculate blue (positive) density map blobs. Default option if red/green not selected.  However, this option uses a LOT OF MEMORY at 1.5sd.
    --radius=<radius>               Radius (in angstroms) around atom or residue to calculate significant discrepancy. [default: 3.5]
    --num-sd=<num-sd>               Number of standard deviation units to use as a significant discrepancy cutoff. Default is 3.0 for atom, residue, green, and red.  Default is 1.5 for blue.
    --type=<type>                   Residue type or atom type to filter by.
    --out-format=<format>           Output file format, available formats: csv, json [default: json].
    --print-validation              Print comparison of median absolute values below 1 sigma for Fo and Fc maps, which should be very similar.
"""

import docopt
import sys
import json
import jsonpickle
import numpy

from . import densityAnalysis
from . import __version__

def main():
    args = docopt.docopt(__doc__, version=__version__)
    if args["--help"]:
        print(__doc__)
        exit(0)

    args["--radius"] = float(args["--radius"])

    if args["--num-sd"] == None:
        args["--num-sd"] = 3.0 if args["--green"] or args["--red"] or args["difference"] else 1.5
    args["--num-sd"] = float(args["--num-sd"])

    if args["--params"]:
        try:
            with open(args["--params"], 'r') as paramsFile:
                params = json.load(paramsFile)
            densityAnalysis.setGlobals(params)
        except:
            sys.exit(str("Error: params file \"") + args["--params"] + "\" does not exist or is not parsable.")

    analyzer = densityAnalysis.fromPDBid(args["<pdbid>"])
    if not analyzer:
        sys.exit("Error: Unable to parse or download PDB entry or associated ccp4 file.")

    jsonType = "json"
    if args["--density"]:
        jsonType = "jsonpickle"
        result = analyzer.densityObj
    elif args["--diff-density"]:
        jsonType = "jsonpickle"
        result = analyzer.diffDensityObj
    elif args["density"]:
        analyzer.aggregateCloud()
        if args["--atom"]:
            headerList = list(map(str, list(analyzer.atomCloudDescriptions.dtype.names) + ['density_electron_ratio']))
            result = [[numpyConverter(element) for element in item] + [analyzer.densityElectronRatio] for item in analyzer.atomCloudDescriptions]
        elif args["--residue"]:
            headerList = densityAnalysis.DensityAnalysis.residueCloudHeader + ['density_electron_ratio']
            result = [list(item) + [analyzer.densityElectronRatio] for item in analyzer.residueCloudDescriptions]
        elif args["--domain"]:
            headerList = densityAnalysis.DensityAnalysis.domainCloudHeader + ['density_electron_ratio']
            result = [list(item) + [analyzer.densityElectronRatio] for item in analyzer.domainCloudDescriptions]
    elif args["difference"]:
        if args["--atom"]:
            headerList = densityAnalysis.DensityAnalysis.atomRegionDiscrepancyHeader
            result = analyzer.calculateAtomRegionDiscrepancies(args["--radius"], args["--num-sd"], args["--type"])
        elif args["--residue"]:
            headerList = densityAnalysis.DensityAnalysis.residueRegionDiscrepancyHeader
            result = analyzer.calculateResidueRegionDiscrepancies(args["--radius"], args["--num-sd"], args["--type"])
        elif args["--symmetry-atom"]:
            headerList = densityAnalysis.DensityAnalysis.symmetryAtomRegionDiscrepancyHeader
            result = analyzer.calculateSymmetryAtomRegionDiscrepancies(args["--radius"], args["--num-sd"], args["--type"])
            for atomInfo in result:
                atomInfo[4] = [val for val in atomInfo[4]]
                atomInfo[5] = [float(val) for val in atomInfo[5]]
    elif args['blob']:
        headerList = densityAnalysis.DensityAnalysis.blobStatisticsHeader
        result = []
        if args["--green"]:
            greenBlobList = analyzer.diffDensityObj.createFullBlobList(analyzer.diffDensityObj.meanDensity + args["--num-sd"] * analyzer.diffDensityObj.stdDensity)
            result.extend(analyzer.calculateAtomSpecificBlobStatistics(greenBlobList))
        if args["--red"]:
            redBlobList = analyzer.diffDensityObj.createFullBlobList(-1 * (analyzer.diffDensityObj.meanDensity + args["--num-sd"] * analyzer.diffDensityObj.stdDensity))
            result.extend(analyzer.calculateAtomSpecificBlobStatistics(redBlobList))
        if not args["--green"] and not args["--red"]: # args["--blue"] by default
            blueBlobList = analyzer.densityObj.createFullBlobList(analyzer.densityObj.meanDensity + args["--num-sd"] * analyzer.densityObj.stdDensity)
            result.extend(analyzer.calculateAtomSpecificBlobStatistics(blueBlobList))
        for blobInfo in result:
            blobInfo[9]  = [ val for val in blobInfo[9] ]
            blobInfo[10] = [ float(val) for val in blobInfo[10] ]
            blobInfo[11] = [ float(val) for val in blobInfo[11] ]
    elif args["statistics"]:
        if args["--print-validation"]:
            (medianAbsFo,medianAbsFc) = analyzer.medianAbsFoFc()
            print("Median abs Fo(<1sd):", medianAbsFo, "Median abs Fc(<1sd):", medianAbsFc,"Relative Difference:",(medianAbsFo - medianAbsFc)/max(medianAbsFo,medianAbsFc))

        if args["--residue"]:
            headerList = analyzer.residueMetricsHeaderList
            result = analyzer.residueMetrics()
        elif args["--atom"]:
            headerList = analyzer.atomMetricsHeaderList
            result = analyzer.atomMetrics()
            for atomInfo in result:
                atomInfo[4] = [x for x in atomInfo[4]]
                atomInfo[5] = [float(x) for x in atomInfo[5]]

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

def numpyConverter(obj):
    """Converts numpy objects to standard Python types, otherwise just return the object.

    :param obj:
    :return: object
    """
    if isinstance(obj, numpy.integer):
        return int(obj)
    elif isinstance(obj, numpy.floating):
        return float(obj)
    elif isinstance(obj, numpy.ndarray):
        return [numpyConverter(item) for item in obj]
    else:
        return obj
