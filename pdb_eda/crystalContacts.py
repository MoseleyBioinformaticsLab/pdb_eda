#!/usr/bin/python3
"""
pdb_eda (crystal) contacts analysis mode command-line interface
	Analyzes atoms in a PDB entry for the closest crystal contacts.
	This mode requires the pymol package and associated python library to be installed.

Usage:
	pdb_eda contacts -h | --help
	pdb_eda contacts <pdbid> <out-file> [--distance=<cutoff>] [--symmetry-atoms] [--include-pdbid] [--out-format=<format>]

Options:
    -h, --help             Show this screen.
    <pdbid>                The PDB ID to download and analyze.
    <out-file>             Output filename. "-" will write to standard output.
    --distance=<cutoff>    Distance cutoff in angstroms for detecting crystal contacts [default: 5.0].
    --symmetry-atoms       Calculate crystal contacts to symmetry atoms too.
    --include-pdbid        Include PDB ID at the beginning of each result.
    --out-format=<format>  Output file format, available formats: csv, json [default: json].
"""

import scipy.spatial.distance
import numpy as np
import sys
import os
import docopt
import json

from . import densityAnalysis
from . import __version__

try:
	import pymol
	pymolImported = True
except ImportError:
	pymolImported = False

def main():
	if not pymolImported:
		print("Error: pymol module did not load.")
		print(__doc__)
		exit(0)

	args = docopt.docopt(__doc__, version=__version__)

	if args["--help"]:
		print(__doc__)
		exit(0)

	args["--distance"] = float(args["--distance"])
	args["<pdbid>"] = args["<pdbid>"].lower()

	analyzer = densityAnalysis.fromPDBid(args["<pdbid>"], mmcif=True)
	if not analyzer:
		raise RuntimeError("Error: Unable to parse or download PDB entry or associated ccp4 file.")

	mmcif_file = densityAnalysis.pdbfolder + args["<pdbid>"] + '.cif.gz'

	crystalNeighborCoords = simulateCrystalNeighborCoordinates(mmcif_file, args["--distance"])

	if args["--symmetry-atoms"]:
		atoms = analyzer.symmetryAtoms
		contacts = findCoordContacts(analyzer.symmetryAtomCoords, crystalNeighborCoords, args["--distance"])
	else:
		atoms = list(analyzer.biopdbObj.get_atoms())
		contacts = findCoordContacts(np.asarray([atom.coord for atom in atoms]), crystalNeighborCoords, args["--distance"])

	headerList = ['model', 'chain', 'residue_number', 'residue_name', "atom_name", "occupancy", "symmetry", "xyz", "crystal_contact_distance"]
	result = []
	for index,contactDistance in contacts:
		atom = atoms[index]
		result.append([atom.parent.parent.parent.id, atom.parent.parent.id, atom.parent.id[1], atom.parent.resname, atom.name, atom.get_occupancy(),
					   [x for x in atom.symmetry] if args["--symmetry-atoms"] else [0,0,0,0], [float(c) for c in atom.coord], contactDistance])

	if args["--include-pdbid"]:
		headerList = ["pdbid"] + headerList
		result = [[analyzer.pdbid] + row for row in result]

	with open(args["<out-file>"], 'w') if args["<out-file>"] != "-" else sys.stdout as outFile:
		if args["--out-format"] == 'csv':
			csvResults = [','.join(map(str, row)) for row in [headerList] + result]
			print(*csvResults, sep='\n', file=outFile)
		else:
			jsonResults = [dict(zip(headerList, row)) for row in result]
			print(json.dumps(jsonResults, indent=2, sort_keys=True), file=outFile)


def findCoordContacts(coordList1, coordList2, distanceCutoff=5.0):
	"""Find contacts in coordList1 to coordList2 at the given distance cutoff.

	:param coordList1: list of coordinates.
	:type coordList1: :py:class:`list`
	:param coordList2: list of coordinates.
	:type coordList2: py:class:`list`
	:param distanceCutoff: distance cutoff., defaults to 5.0
	:type distanceCutoff: :py:class:`float`

	:return: contactList of index,minDistance tuples.
	:rtype: :py:class:`list`
	"""
	distances = scipy.spatial.distance.cdist(coordList1, coordList2)
	return [(index,minDistance) for index,minDistance in enumerate(np.min(distances[x]) for x in range(len(coordList1))) if minDistance <= distanceCutoff]


def simulateCrystalNeighborCoordinates(filename, distanceCutoff=5.0):
	"""RETURN a list of atom coordinates of the simulated crystal environment surrounding the X-Ray
	Diffraction asymmetric unit (excluding heteroatoms). Requires a file path instead of a
	structure because the bulk of this is handled by Pymol.
	NOTE: This will only work with PDB structures resolved with X-RAY DIFFRACTION.

	:param filename:
	:type filename: :py:class:`str`
	:param distanceCutoff: distance cutoff., defaults to 5.0
	:type distanceCutoff: :py:class:`float`

	:return: coordList
	:rtype: :py:class:`list`
	"""
	# Launch pymol.
	pymol.pymol_argv = ['pymol', '-qc']
	pymol.finish_launching()

	# Set absolute path and name of PDB entry.
	spath = os.path.abspath(filename)
	sname = spath.split('/')[-1].split('.')[0]
	asym_unit = "asym_unit"

	# Load Structure.
	pymol.cmd.load(spath)
	pymol.cmd.disable("all")
	pymol.cmd.enable(sname)
	pymol.cmd.create(asym_unit, "polymer")

	# Generate local crystal environment using symexp; delete original structures.
	pymol.cmd.symexp("neighbor", asym_unit, asym_unit, distanceCutoff)
	pymol.cmd.delete(sname)
	pymol.cmd.delete(asym_unit)

	# Return atom coordinates of the neighboring structures.
	myspace = {"coordinates": []}
	pymol.cmd.iterate_state(1, "all", "coordinates.append([x,y,z])", space=myspace)

	pymol.cmd.reinitialize()
	return myspace["coordinates"]

