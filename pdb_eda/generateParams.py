"""
generateParams.py

Usage:
    pdb_eda generate residue-list




ftp://ftp.wwpdb.org/pub/pdb/data/monomers/components.cif
http://remediation.wwpdb.org/data/ccd

"""

import CifFile
import docopt
import urllib.request
import os
#from . import __version__

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
    #args = docopt.docopt(__doc__, version=__version__)
    args = docopt.docopt(__doc__)
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
        if "_chem_comp_atom.atom_id" in residue and "_chem_comp_atom.charge" in residue and "_chem_comp_atom.type_symbol" in residue and "_chem_comp_atom.pdbx_leaving_atom_flag" in residue and "_chem_comp_atom.pdbx_aromatic_flag" in residue:
            atoms = { name : {"name": name, "charge": charge, "element" : element, "leaving" : leaving, "aromatic" : aromatic, "bonds" : [] } for (name, charge, element, leaving, aromatic) in
                      zip(residue["_chem_comp_atom.atom_id"], residue["_chem_comp_atom.charge"], residue["_chem_comp_atom.type_symbol"], residue["_chem_comp_atom.pdbx_leaving_atom_flag"], residue["_chem_comp_atom.pdbx_aromatic_flag"]) }
            if "_chem_comp_bond.atom_id_1" in residue and "_chem_comp_bond.atom_id_2" in residue and "_chem_comp_bond.value_order" in residue and
                    "_chem_comp_bond.pdbx_aromatic_flag" in residue and "_chem_comp_bond.pdbx_stereo_config" in residue:
                for (atomName1, atomName2, bondType,isAromatic,stereo) in
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