#!/usr/bin/python3
"""
pdb_eda command-line interface
    The PDB Electron Density Analysis package has facilities for analyzing PDB entries and their associated 2Fo-Fc and Fo-Fc electron density maps.

Usage:
    pdb_eda -h | --help     for this screen.
    pdb_eda --full-help     help documentation on all modes.
    pdb_eda --version       for the version of pdb_eda.
    pdb_eda single ...      for single structure analysis mode. (Most useful command line mode).
    pdb_eda multiple ...    for multiple structure analysis mode. (Second most useful command line mode).
    pdb_eda contacts ...    for crystal contacts analysis mode.  (Third most useful command line mode).
    pdb_eda generate ...    for generating starting parameters file that then needs to be optimized. (Rarely used mode).
    pdb_eda optimize ...    for parameter optimization mode. (Rarely used mode).

For help on a specific mode, use the mode option -h or --help.
For example:
    pdb_eda single --help   for help documentation about single structure analysis mode.
"""

import sys
from . import singleStructure
from . import multipleStructures
from . import crystalContacts
from . import generateParams
from . import optimizeParams
from . import __version__

def main():
    if len(sys.argv) > 1 and sys.argv[1] == "single":
        singleStructure.main()
    elif len(sys.argv) > 1 and sys.argv[1] == "multiple":
        multipleStructures.main()
    elif len(sys.argv) > 1 and sys.argv[1] == "optimize":
        optimizeParams.main()
    elif len(sys.argv) > 1 and sys.argv[1] == "contacts":
        crystalContacts.main()
    elif len(sys.argv) > 1 and sys.argv[1] == "generate":
        generateParams.main()
    elif len(sys.argv) > 1 and (sys.argv[1] == "--version" or sys.argv[1] == "-v") :
        print("Version: ",__version__)
    elif len(sys.argv) > 1 and sys.argv[1] == "--full-help":
        print(__doc__)
        print("-"*80)
        print(singleStructure.__doc__)
        print("-"*80)
        print(multipleStructures.__doc__)
        print("-"*80)
        print(crystalContacts.__doc__)
        print("-"*80)
        print(generateParams.__doc__)
        print("-"*80)
        print(optimizeParams.__doc__)
    else:
        print(__doc__)


if __name__ == "__main__": # this is hidden from the console script created by pip, since main() is called directly.
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