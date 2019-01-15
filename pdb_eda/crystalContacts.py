#!/usr/bin/python3

"""crystal_contacts.py
Python3 script for identifying 'crystal contact' atoms in a PDB entry.  If run as a script, saves a jsonpickled list of Biopython Atom objects.

If imported as a module, simply use the get_contact_atoms(mmcif_file) function for any PDB entry.

Usage:
	crystal_contacts.py <pdb_mmcif_file>
	crystal_contacts.py -h

Arguments:
	<pdb_mmcif> - /path/to/pdb/1xyz.cif
"""

import Bio.PDB.MMCIFParser
import itertools
import scipy.spatial.distance
import numpy as np
import pymol as pm
import os
import gzip
import jsonpickle
import docopt


def main(args):
	mmcif_file = args["<pdb_mmcif_file>"]

	entry_id = mmcif_file.split('/')[-1].split('.')[0]

	# Get list of contact atoms as Biopython Atom objects.
	contact_atoms = get_contact_atoms(mmcif_file)

	# Save list of atoms as a json pickle.
	json_atom_list = jsonpickle.encode(contact_atoms)
	with open(entry_id + "_contact_atoms.json", "w") as json_file:
		json_file.write(json_atom_list)

	exit(0)


def _pymol_test(mmcif_file):
	"""For verifying that this method of identifying crystal contacts is reliable."""

	contact_atoms = get_contact_atoms(mmcif_file)

	pm.cmd.reinitialize()
	pm.cmd.load(mmcif_file)
	pm.cmd.show_as("surface")
	pm.cmd.color("orange", "all")

	for atom in contact_atoms:
		pm.cmd.color("red", "///{0}/{1}/{2}".format(atom.get_full_id()[2], atom.get_full_id()[3][1], atom.get_full_id()[4][0]))


def get_contact_atoms(mmcif_file):
	"""Return a list of Biopython Atom objects which are involved in crystal contacts for the given PDB entry (mmcif format)."""

	parser = Bio.PDB.MMCIFParser()

	# Open with gzip if compressed.
	if mmcif_file.endswith("gz"):
		mmcif_filehandle = gzip.open(mmcif_file, "rt")
	else:
		mmcif_filehandle = open(mmcif_file, "rt")

	# Parse mmcif file into Biopython Structure.
	structure = parser.get_structure("entry", mmcif_filehandle)

	# Get atoms of asymmetric unit (exluding heteroatoms).
	asym_unit_atoms = [atom for atom in structure.get_atoms() if atom.get_full_id()[3][0] == " "]
	asym_unit_coords = [atom.get_coord() for atom in asym_unit_atoms]

	# Get atoms of neighboring crystal environment.
	crystal_neighbor_coords = simulate_crystal_coordinates(mmcif_file)

	# Determine which atoms are 'in contact' with the neighboring crystal units.
	contacts = find_atomic_contacts([asym_unit_coords, crystal_neighbor_coords])
	contact_atoms = [atom for atom, is_contact in zip(asym_unit_atoms, contacts[0]) if is_contact]

	return contact_atoms


def find_atomic_contacts(macromolecules, contact_cutoff=5):
	"""Given a list of atomic coordinate lists, return a "parallel" list of boolean atomic contact lists.
	Parameters:
		macromolecules - a list of lists; each internal list contains all atoms in a single molecular entity.
		contact_cutoff - atoms within this distance (in angstroms) of another molecule are considered 'atomic contacts'."""

	print("Running cdist...")

	# Build an array of atomic distance matrices for all pairwise combinations (i, j) of macromolecules in macromolecules_coordinates.
	dist_matrices = np.array(
		[scipy.spatial.distance.cdist(mol_i, mol_j) for mol_i, mol_j in itertools.combinations(macromolecules, 2)])

	print("Done running cdist.")

	# Build an array of intermolecular atomic contact maps for each pair of macromolecules.
	pairwise_contact_maps = np.array(
		[[np.less(row, contact_cutoff) for row in dist_matrix] for dist_matrix in dist_matrices])

	# Initially all atoms are assumed to be non-contact atoms (False), and updated if indicated by a contact map in the below iteration.
	contact_axes = [[False for atom in mol] for mol in macromolecules]

	# Keep track of which macromolecules are represented in each combinatorial distance matrix by their index in 'macromolecules'.
	index_pairs = [index_pair for index_pair in itertools.combinations(range(0, len(macromolecules)), 2)]

	# Identify which atoms make at least one intermolecule contact for each macromolecule in each contact map.
	for cmap, index_pair in zip(pairwise_contact_maps, index_pairs):
		# Find truth values for atomic contacts for each macromolecule in the current contact map.
		axis_1_contact = np.squeeze(np.asarray(cmap.any(1)))
		axis_0_contact = np.squeeze(np.asarray(cmap.any(0)))

		# Update truth values for atomic contacts for each macromolecule.  If atom contact is True for ANY contact map, then the atom is a contact atom.
		contact_axes[index_pair[0]] = [True if (was_true or now_true) else False for was_true, now_true in
									   zip(contact_axes[index_pair[0]], axis_1_contact)]
		contact_axes[index_pair[1]] = [True if (was_true or now_true) else False for was_true, now_true in
									   zip(contact_axes[index_pair[1]], axis_0_contact)]

	return contact_axes


def simulate_crystal_coordinates(filename):
	"""RETURN a list of atom coordinates of the simulated crystal environment surrounding the X-Ray
	Diffraction asymmetric unit (excluding heteroatoms). Requires a file path instead of a 
	structure because the bulk of this is handled by Pymol.
	NOTE: This will only work with PDB structures resolved with X-RAY DIFFRACTION.
	Parameters:
		filename - /path/to/pdb/file"""

	file_path = filename

	# Launch pymol.
	pm.pymol_argv = ['pymol', '-qc']
	pm.finish_launching()

	# Set absolute path and name of PDB entry.
	spath = os.path.abspath(file_path)
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


if __name__ == "__main__":
	args = docopt.docopt(__doc__)
	main(args)
