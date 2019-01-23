# !/usr/bin/python3

"""
pdb_eda command-line interface

Usage:
    pdb_eda -h | --help
    pdb_eda --version
    pdb_eda (<pdbid> <out-file>) [--radii-param=<paramfile>] [--atom] [--residue] [--chain] [--out-format=<format>]
    pdb_eda (<pdbid> <out-file>) [--radii-param=<paramfile>] [--green | --red | --all] [--stats] [--out-format=<format>]
    pdb_eda (<pdbid> <out-file>) [--radii-param=<paramfile>] [--symmetry-atoms]

Options:
    -h, --help                      Show this screen.
    --version                       Show version.
    <pdbid>                         The PDB id
    <out-file>                      Output file name
    --radii-param=<paramfile>       Radii parameters. [default: ../conf/optimized_radii_slope_param.json]
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

import jsonpickle
from . import densityAnalysis

## Final set of radii derived from optimization
radiiDefault = {'C_single': 0.84, 'C_double': 0.67, 'C_intermediate': 0.72, 'C_single_bb': 0.72, 'C_double_bb': 0.61, 
                'O_single': 0.80, 'O_double': 0.77, 'O_intermediate': 0.82, 'O_double_bb': 0.71, 
                'N_single': 0.95, 'N_intermediate': 0.77, 'N_single_bb': 0.7, 
                'S_single': 0.75}

slopesDefault = {'C_single': -0.33373359301010996, 'C_double': -0.79742538209281033, 'C_intermediate': -0.46605044936397311, 'C_single_bb': -0.44626029983492604, 'C_double_bb': -0.53313491315106321, 
                'O_single': -0.61612691527685515, 'O_double': -0.69081386892018048, 'O_intermediate': -0.64900152439990022, 'O_double_bb': -0.59780395867126568, 
                'N_single': -0.50069328952200343, 'N_intermediate': -0.5791543104296577, 'N_single_bb': -0.4905830761073553, 
                'S_single': -0.82192528851026947}


def fromPDBid(args):
    pdbid = args["<pdbid>"]
    filename = args["<out-file>"]
    #if not args["<radii-param>"]:
    #    pass

    analyser = densityAnalysis.fromPDBid(pdbid)
    result = []
    if args["--atom"] or args["--residue"] or args["--chain"]:
        analyser.aggregateCloud(radiiDefault, slopesDefault, atomL=True, residueL=True, chainL=True)
        if args["--atom"]:
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
            result = analyser.calcAtomBlobDists()
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
