# !/usr/bin/python3

import os
import json
import jsonpickle
from . import densityAnalysis


def main(args):
    pdbid = args["<pdbid>"]
    filename = args["<out-file>"]

    paramsPath = os.path.join(os.path.dirname(__file__), args["--radii-param"])
    with open(paramsPath, 'r') as fh:
        params = json.load(fh)
    radii = params['radii']
    slopes = params['slopes']

    analyser = densityAnalysis.fromPDBid(pdbid)
    result = []
    if args["--density-map"]:
        result = analyser.densityObj
    elif args["--diff-density-map"]:
        result = analyser.diffDensityObj


    if args["--atom"] or args["--residue"] or args["--chain"]:
        analyser.aggregateCloud(radii, slopes, atomL=True, residueL=True, chainL=True)
        if args["--atom"]:
            result.append(','.join(map(str, list(analyser.atomList) + ['chainMedian'])))
            for item in analyser.atomList.values.tolist():
                result.append(','.join(map(str, item + [analyser.chainMedian])))
        if args["--residue"]:
            for item in analyser.residueList:
                result.append(','.join(map(str, item + [analyser.chainMedian])))
        if args["--chain"]:
            for item in analyser.chainList:
                result.append(','.join(map(str, item + [analyser.chainMedian])))

    if args["--green"] or args["--red"] or args["--all"]:
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

    if args["--symmetry-atoms"]:
        analyser.calcSymmetryAtoms()
        result = analyser.symmetryAtoms


    with open(filename, 'w') as fh:
        if args["--out-format"] == 'json':
            result = jsonpickle.encode(result)
            fh.write(result)
        else:
            print(*result, sep='\n', file=fh)
