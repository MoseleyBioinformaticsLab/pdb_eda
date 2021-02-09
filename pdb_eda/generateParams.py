"""
generateParams.py

Usage:
    pdb_eda optimize -h | --help
    pdb_eda generate color [--params=<params-file>]
    pdb_eda generate residue-report <pdbid-file> <out-jsonfile>

Options:
    --params=<params-file>      Overriding parameters file that includes radii, slopes, etc. [default: ]



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
defaultParamsFilepath = os.path.join(os.path.dirname(__file__), 'conf/optimized_params.json')

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

    paramsFilepath = args["--params"] if args["--params"] else defaultParamsFilepath
    try:
        with open(paramsFilepath, 'r') as paramsFile:
            params = json.load(paramsFile)
    except:
        sys.exit(str("Error: params file \"") + paramsFilepath + "\" does not exist or is not parsable.")


    if args["color"]:
        for residue in componentsInfo["residues"].values():
            residueName = residue["name"].upper()
            if residueName in params["residue_electrons"]:
                residue["estimated_electrons"] = 0
                for atom in residue["atoms"].values():
                    atom["num_bound_hydrogens"] = sum(1 for atomName,bondType,aromatic,stereo in atom["bonds"] if atomName in residue["atoms"] and residue["atoms"][atomName]["element"] == "H" and residue["atoms"][atomName]["leaving"] == atom["leaving"])
                    if atom["element"] in params["element_electrons"]:
                        atom["estimated_electrons"] = params["element_electrons"][atom["element"]] + atom["num_bound_hydrogens"] - float(atom["charge"])
                        residue["estimated_electrons"] += atom["estimated_electrons"] if atom["leaving"] != "Y" else 0
                    atom["color"] = atom["name"] + "." + atom["element"] + "." + atom["aromatic"]
                    atom["element_color"] = atom["element"] + "." + atom["aromatic"]
                for atom in residue["atoms"].values():
                    atom["full_bond_colors"] = [ residue["atoms"][atomName]["color"] + "." +  bondType + "." + aromatic for atomName,bondType,aromatic,stereo in atom["bonds"] if atomName in residue["atoms"]]
                    atom["element_bond_colors"] = [ residue["atoms"][atomName]["element_color"] + "." +  bondType + "." + aromatic for atomName,bondType,aromatic,stereo in atom["bonds"] if atomName in residue["atoms"]]
                    atom["full_color"] = atom["color"] + "-" + "_".join(sorted(atom["full_bond_colors"]))
                    atom["full_element_color"] = atom["element_color"] + "#" + "_".join(sorted(atom["element_bond_colors"]))
                print(residueName, "Estimated Electrons:", params["residue_electrons"][residueName], residue["estimated_electrons"])
                for atom in residue["atoms"].values():
                    atomName = residueName + "_" + atom["name"]
                    if atomName in params["full_atom_name_map_electrons"]:
                        print(" ", atomName, "\tEstimated Electrons:", params["full_atom_name_map_electrons"][atomName], atom["estimated_electrons"])
    elif args["residue-report"]:
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



