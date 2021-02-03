#!/usr/bin/python3
"""
pdb_eda (crystal) contacts analysis mode command-line interface
	Analyzes atoms in a PDB entry for the closest crystal contacts.
	This mode requires the pymol package and python library to be installed.

Usage:
	pdb_eda contacts -h | --help
	pdb_eda contacts <pdbid> <out-file> [--distance=<cutoff>] [--symmetry-atoms] [--include-pdbid] [--out-format=<format>]

Arguments:
    <pdbid>				  The PDB ID to download and analyze.
    <out-file>			  Output filename. "-" will write to standard output.
    --distance=<cutoff>	  Distance cutoff in angstroms for detecting crystal contacts. [default: 5.0]
    --symmetry-atoms	  Calculate crystal contacts to symmetry atoms too.
    --include-pdbid		  Include PDB ID at the beginning of each result.
"""

import scipy.spatial.distance
import numpy as np
import pymol as pm
import os
import docopt
import json

from . import densityAnalysis
from . import __version__

def main():
	args = docopt.docopt(__doc__, version=__version__)
	if args["--help"]:
		print(__doc__)
		exit(0)

	args["--distance"] = float(args["--distance"])

	analyzer = densityAnalysis.fromPDBid(args["<pdbid>"], mmcif=True)
	if not analyzer:
		sys.exit("Error: Unable to parse or download PDB entry or associated ccp4 file.")

	mmcif_file = densityAnalysis.pdbfolder + args["<pdbid>"] + '.cif'

	crystalNeighborCoords = simulateCrystalNeighborCoordinates(mmcif_file)

	if args["--symmetry-atoms"]:
		contacts = findCoordContacts(analyzer.symmetryAtomCoords, crystalNeighborCoords, args["--distance"])
		atoms = analyzer.symmetryAtoms
	else:
		contacts = findCoordContacts(analyzer.asymmetryAtomCoords, crystalNeighborCoords, args["--distance"])
		atoms = analyzer.asymmetryAtoms

	headerList = ['chain', 'residue_number', 'residue_name', "atom_name", "min_occupancy", "atom_symmetry", "atom_xyz", "crystal_contact_distance"]
	result = []
	for index,contactDistance in contacts:
		atom = atoms[index]
		result.append([atom.parent.parent.id, atom.parent.id[1], atom.parent.resname, atom.name, atom.get_occupancy(), atom.symmetry, [float (c) for c in atom.coord], contactDistance])

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
	"""
	Find contacts in coordList1 to coordList2 at the given distance cutoff.

	:param :py:obj:`list` coordList1: list of coordinates.
	:param :py:obj:`list` coordList2: list of coordinates.
	:param float distanceCutoff: distance cutoff.
	:return: contactList of index,minDistance tuples.
    :rtype: :py:obj:`list`
	"""
	distances = scipy.spatial.distance.cdist(np.matrix(coordList1), np.matrix(coordList2))
	return [(index,minDistance) for index,minDistance in enumerate(np.min(distances[x]) for x in range(len(coordList1))) if minDistance <= distanceCutoff]


def simulateCrystalNeighborCoordinates(filename):
	"""
	RETURN a list of atom coordinates of the simulated crystal environment surrounding the X-Ray
	Diffraction asymmetric unit (excluding heteroatoms). Requires a file path instead of a
	structure because the bulk of this is handled by Pymol.
	NOTE: This will only work with PDB structures resolved with X-RAY DIFFRACTION.

	:param str filename:
	:return: coordList
    :rtype: :py:obj:`list`
	"""
	# Launch pymol.
	pm.pymol_argv = ['pymol', '-qc']
	pm.finish_launching()

	# Set absolute path and name of PDB entry.
	spath = os.path.abspath(filename)
	sname = spath.split('/')[-1].split('.')[0]
	asym_unit = "asym_unit"

	# Load Structure.
	pm.cmd.load(spath, format="cif")
	pm.cmd.disable("all")
	pm.cmd.enable(sname)
	pm.cmd.create(asym_unit, "polymer")

	# Generate local crystal environment using symexp; delete original structures.
	pm.cmd.symexp("neighbor", asym_unit, asym_unit, 5)
	pm.cmd.delete(sname)
	pm.cmd.delete(asym_unit)

	# Return atom coordinates of the neighboring structures.
	myspace = {"coordinates": []}
	pm.cmd.iterate_state(1, "all", "coordinates.append([x,y,z])", space=myspace)

	pm.cmd.reinitialize()
	return myspace["coordinates"]

