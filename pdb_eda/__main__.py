#!/usr/bin/python3
"""
pdb_eda command-line interface

Usage:
    pdb_eda single ...      for single structure analysis mode.
    pdb_eda multiple ...    for multiple structure analysis mode.
    pdb_eda -h | --help     for this screen.
    pdb_eda --version       for the version of pdb_eda.


For help on a specific mode, use the mode option -h or --help.
For example:
    pdb_eda single --help   for help documentation about single structure analysis mode.
"""

import sys
from . import singleStructure
from . import multipleStructures
from . import __version__


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "single":
        singleStructure.main()
    elif len(sys.argv) > 1 and sys.argv[1] == "multiple":
        multipleStructures.main()
    elif len(sys.argv) > 1 and (sys.argv[1] == "--version" or sys.argv[1] == "-v") :
        print("Version: ",__version__)
    else:
        print(__doc__)


if __name__ == "__main__":
    main()
