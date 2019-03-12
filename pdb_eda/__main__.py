#!/usr/bin/python3

"""
pdb_eda command-line interface

Usage:
    pdb_eda -h | --help
    pdb_eda --version
    pdb_eda single <pdbid> <out-file> [--density-map | --diff-density-map]
    pdb_eda single <pdbid> <out-file> [--radii-param=<paramfile>] [--atom] [--residue] [--chain] [--out-format=<format>]
    pdb_eda single <pdbid> <out-file> [--radii-param=<paramfile>] [--green | --red | --all] [--stats] [--out-format=<format>]
    pdb_eda single <pdbid> <out-file> [--radii-param=<paramfile>] [--symmetry-atoms]
    pdb_eda multiple <pdbid-file> <out-file> [--radii-param=<paramfile>]

Options:
    -h, --help                      Show this screen.
    --version                       Show version.
    single                          Running single-structure mode
    <pdbid>                         The PDB id
    <out-file>                      Output file name
    multiple                        Running multiple-structure mode
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

from . import singleStructure
from . import multipleStructures
from . import __version__


def main(args):
    if args["single"]:
        singleStructure.main(args)
    elif args["multiple"]:
        multipleStructures.main(args)


if __name__ == "__main__":
    args = docopt.docopt(__doc__, version=__version__)
    main(args)
