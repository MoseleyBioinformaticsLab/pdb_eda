"""
generateParams.py

Usage:
    pdb_eda optimize -h | --help
    pdb_eda generate residue-report <pdbid-file> <out-jsonfile>




ftp://ftp.wwpdb.org/pub/pdb/data/monomers/components.cif
http://remediation.wwpdb.org/data/ccd

"""

import CifFile
import docopt
import urllib.request
import os
import collections
import json

from . import __version__
from . import densityAnalysis

componentsFilename = "components.cif"

elementElectrons: {
    "AL": 13,
    "AU": 79,
    "BA": 56,
    "BR": 35,
    "C": 6,
    "CA": 20,
    "CD": 48,
    "CL": 17,
    "CO": 27,
    "CR": 24,
    "CS": 55,
    "CU": 29,
    "F": 9,
    "FE": 26,
    "HG": 80,
    "I": 53,
    "K": 19,
    "LI": 3,
    "MG": 12,
    "MN": 25,
    "MO": 42,
    "N": 7,
    "NA": 11,
    "NI": 28,
    "O": 8,
    "P": 15,
    "PT": 78,
    "RU": 44,
    "S": 16,
    "SR": 38,
    "V": 23,
    "W": 74,
    "Y": 39,
    "YB": 70,
    "ZN": 30
  }


def main():
    args = docopt.docopt(__doc__, version=__version__)
    if args["--help"]:
        print(__doc__)
        exit(0)

    if os.path.isfile("components_info.json"):
        with open("components_info.json", 'r') as jsonFile:
            componentsInfo = json.load(jsonFile)
    else:
        componentsInfo = processComponents()

    if not os.path.isfile("components_info.json"):
        with open("components_info.json", 'w') as outFile:
            print(json.dumps(componentsInfo, indent=2, sort_keys=True), file=outFile)

    if args["residue-report"]:
        try:
            pdbids = []
            with open(args["<pdbid-file>"], "r") as textFile:
                for pdbid in textFile:
                    pdbids.append(pdbid[0:4])
        except:
            sys.exit(str("Error: PDB IDs file \"") + args["<pdbid-file>"] + "\" does not exist or is not parsable.")

        totalResidueAtoms = collections.defaultdict(int)
        totalElements = collections.defaultdict(int)
        totalResidues = collections.defaultdict(int)
        pdbidInfo = collections.defaultdict(dict)
        for pdbid in pdbids:
            analyzer = densityAnalysis.fromPDBid(pdbid)
            if not analyzer:
                continue

            pdbidInfo[pdbid]["properties"] = { property : value for (property,value) in analyzer.biopdbObj.header.items() }
            pdbidInfo[pdbid]["properties"]["resolution"] = analyzer.pdbObj.header.resolution
            pdbidInfo[pdbid]["properties"]['voxel_volume'] = analyzer.densityObj.header.unitVolume
            pdbidInfo[pdbid]["properties"]['space_group'] = analyzer.pdbObj.header.spaceGroup

            pdbidInfo[pdbid]["residue_atom_counts"] = collections.Counter(densityAnalysis.residueAtomName(atom) for residue in analyzer.biopdbObj.get_residues() for atom in residue.get_atoms())
            for name,count in pdbidInfo[pdbid]["residue_atom_counts"].items():
                totalResidueAtoms[name] += count

            pdbidInfo[pdbid]["element_counts"] = collections.Counter(atom.element for residue in analyzer.biopdbObj.get_residues() for atom in residue.get_atoms())
            for name,count in pdbidInfo[pdbid]["element_counts"].items():
                totalElements[name] += count

            pdbidInfo[pdbid]["residue_counts"] = collections.Counter(residue.resname for residue in analyzer.biopdbObj.get_residues())
            for name,count in pdbidInfo[pdbid]["residue_counts"].items():
                totalResidues[name] += count

        with open(args["<out-jsonfile>"], 'w') if args["<out-jsonfile>"] != "-" else sys.stdout as outFile:
            jsonOutput = { "pdbid" : pdbidInfo, "residue_atom_counts" : totalResidueAtoms, "residue_counts" : totalResidues, "element_counts" : totalElements }
            print(json.dumps(jsonOutput, indent=2, sort_keys=True), file=outFile)



def processComponents():
    if not os.path.isfile(componentsFilename):
        urllib.request.urlretrieve("ftp://ftp.wwpdb.org/pub/pdb/data/monomers/components.cif", componentsFilename)

    print("Parsing the components.cif file.  This takes a few minutes.")
    components = CifFile.ReadCif(componentsFilename)

    print("Processing the components.")
    residues = {}
    errors = set()
    for residueName in components.keys():
        residue = components[residueName]
        if "_chem_comp_atom.atom_id" in residue and "_chem_comp_atom.charge" in residue and "_chem_comp_atom.type_symbol" in residue and "_chem_comp_atom.pdbx_leaving_atom_flag" in residue and \
                "_chem_comp_atom.pdbx_aromatic_flag" in residue:
            atoms = { name : {"name": name, "charge": charge, "element" : element, "leaving" : leaving, "aromatic" : aromatic, "bonds" : [] } for (name, charge, element, leaving, aromatic) in \
                      zip(residue["_chem_comp_atom.atom_id"], residue["_chem_comp_atom.charge"], residue["_chem_comp_atom.type_symbol"], residue["_chem_comp_atom.pdbx_leaving_atom_flag"], residue["_chem_comp_atom.pdbx_aromatic_flag"]) }
            if "_chem_comp_bond.atom_id_1" in residue and "_chem_comp_bond.atom_id_2" in residue and "_chem_comp_bond.value_order" in residue and \
                "_chem_comp_bond.pdbx_aromatic_flag" in residue and "_chem_comp_bond.pdbx_stereo_config" in residue:
                for (atomName1, atomName2, bondType,isAromatic,stereo) in \
                    zip(residue["_chem_comp_bond.atom_id_1"],residue["_chem_comp_bond.atom_id_2"],residue["_chem_comp_bond.value_order"], residue["_chem_comp_bond.pdbx_aromatic_flag"],residue["_chem_comp_bond.pdbx_stereo_config"]):
                    if atomName1 in atoms:
                        atoms[atomName1]["bonds"].append((atomName2,bondType,isAromatic,stereo))
                    else:
                        errors.add(residueName)
                    if atomName2 in atoms:
                        atoms[atomName2]["bonds"].append((atomName1,bondType,isAromatic,stereo))
                    else:
                        errors.add(residueName)
            residues[residueName] = { "name" : residueName, "atoms" : atoms }

    return { "residues" : residues, "errors" : list(errors) }




if __name__ == "__main__":
    main()