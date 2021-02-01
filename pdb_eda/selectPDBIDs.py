#!usr/bin/python3
"""
pdb_eda select PDB IDs mode command-line interface
    Select pdbids (for parameter optimization) based on selection criteria.

    This is used after a set of PDB entries have been analyzed in multiple structure analysis mode and before the optimize mode.
    Also, it is good to use the --time-out option in the multiple structure analysis mode to select entries that analyze quickly.

Usage:
    pdb_eda select -h | --help
    pdb_eda select <pdb-results-file> <out-file> [options]

Options:
    -h, --help                          Show this screen.
    <out-file>                          Output filename. "-" will write to standard output.
    <pdb-results-file>                  JSON File name that contains PDB entry results generated from multiple mode. "-" will read from standard input.
    --out-format=<format>               Output file format, available formats: txt, json [default: txt].
    --max-resolution=<max-resolution>   Maximum x-ray crystallographic resolution to allow. [default: 3.5]
    --min-resolution=<min-resolution>   Minimum x-ray crystallographic resolution to allow. [default: 0]
    --min-atoms=<min-atoms>             Minimum number of atoms analyzed (numAtomsAnalyzed) required. [default: 1000]
    --min-residues=<min-residues>       Minimum number of residues required of the given residue types separated by commas. [default: ]
    --min-elements=<min-elements>       Minimum number of atoms required of the given elements separated by commas. [default: ]

Minimum Residues Example:
   --min-residues=20
      Requires at least 20 residues
   --min-residues=10,ALA,ARG,ASN,ASP,CYS,GLN,GLU,GLY,HIS,ILE,LEU,LYS,MET,PHE,PRO,SER,THR,TRP,TYR,VAL
      Requires at least 10 of the common amino acids.

Minimum Elements Example:
   --min-elements=20
      Requires at least 20 atoms of any given element.
   --min-elements=2,MG
      Requires at least 2 magnesium atoms.
"""
import json
import docopt
import sys

from . import __version__


def main():
    args = docopt.docopt(__doc__, version=__version__)
    if args["--help"]:
        print(__doc__)
        exit(0)

    args["--max-resolution"] = float(args["--max-resolution"])
    args["--min-resolution"] = float(args["--min-resolution"])
    args["--min-atoms"] = int(args["--min-atoms"])

    if "," in args["--min-residues"]:
        allowedResidueTypes = args["--min-residues"].split(",")
        args["--min-residues"] = float(allowedResidueTypes.pop(0))
        allowedResidueTypes = set(allowedResidueTypes)
    else:
        args["--min-residues"] = float(args["--min-residues"])
        allowedResidueTypes = set()

    if "," in args["--min-elements"]:
        allowedElements = args["--min-elements"].split(",")
        args["--min-elements"] = float(allowedElements.pop(0))
        allowedElements = set(allowedElements)
    else:
        args["--min-elements"] = float(args["--min-elements"])
        allowedElements = set()

    try:
        with open(args["<pdb-results-file>"], "r") if args["<pdb-results-file>"] != "-" else sys.stdin as jsonFile:
            pdbResults = json.load(jsonFile)
    except:
        sys.exit(str("Error: PDB results file \"") + args["<pdb-results-file>"] + "\" does not exist or is not parsable.")

    pdbids = [ pdbid for pdbid in pdbResults if
               pdbResults[pdbid]["stats"]["num_atoms_analyzed"] >= args["--min-atoms"] and
               float(pdbResults[pdbid]["stats"]["resolution"]) >= args["--min-resolution"] and
               float(pdbResults[pdbid]["stats"]["resolution"]) <= args["--max-resolution"] and
               sum(count for (residueType,count) in pdbResults[pdbid]["properties"]["residue_counts"].items() if not allowedResidueTypes or residueType in allowedResidueTypes ) >= args["--min-residues"] and
               sum(count for (element, count) in pdbResults[pdbid]["properties"]["element_counts"].items() if not allowedElements or element in allowedElements) >= args["--min-elements"]
             ]

    if args["--out-format"] == "json":
        with open(args['<out-file>'], "w") if args["<out-file>"] != "-" else sys.stdout as jsonFile:
            print(json.dumps(pdbids, indent=2, sort_keys=True), file=jsonFile)
    else:
        with open(args['<out-file>'], "w") if args["<out-file>"] != "-" else sys.stdout as txtFile:
            print("\n".join(pdbids), file=txtFile)

