#!/usr/bin/python3
"""
pdb_eda command-line interface

Usage:
    pdb_eda -h | --help     for this screen.
    pdb_eda --full-help     help documentation on all modes.
    pdb_eda --version       for the version of pdb_eda.
    pdb_eda single ...      for single structure analysis mode.
    pdb_eda multiple ...    for multiple structure analysis mode.
    pdb_eda optimize ...    for radii and slope parameter optimization mode.
    pdb_eda select ...      for PDB entry selection mode.

For help on a specific mode, use the mode option -h or --help.
For example:
    pdb_eda single --help   for help documentation about single structure analysis mode.
"""

import sys
from . import singleStructure
from . import multipleStructures
from . import radiiOptimize
from . import selectPDBIDs
from . import __version__

def main():
    if len(sys.argv) > 1 and sys.argv[1] == "single":
        singleStructure.main()
    elif len(sys.argv) > 1 and sys.argv[1] == "multiple":
        multipleStructures.main()
    elif len(sys.argv) > 1 and sys.argv[1] == "optimize":
        radiiOptimize.main()
    elif len(sys.argv) > 1 and sys.argv[1] == "select":
        selectPDBIDs.main()
    elif len(sys.argv) > 1 and (sys.argv[1] == "--version" or sys.argv[1] == "-v") :
        print("Version: ",__version__)
    elif len(sys.argv) > 1 and sys.argv[1] == "--full-help":
        print(__doc__)
        print("-"*80)
        print(singleStructure.__doc__)
        print("-"*80)
        print(multipleStructures.__doc__)
        print("-"*80)
        print(radiiOptimize.__doc__)
        print("-"*80)
        print(selectPDBIDs.__doc__)
    else:
        print(__doc__)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--profile":
        sys.argv.pop(1)
        import cProfile
        profiler = cProfile.Profile()
        profiler.enable()
        main()
        profiler.disable()
        profiler.print_stats()
    else:
        main()