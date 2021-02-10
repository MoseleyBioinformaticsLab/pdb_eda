"""
generateParams.py

Usage:
    pdb_eda optimize -h | --help
    pdb_eda generate testing
    pdb_eda generate atom-type <out-jsonfile> [--residues=<comma-separated-residues>] [--allow-errors]
    pdb_eda generate prevalence <pdbid-file> <out-jsonfile>
    pdb_eda generate parameters <in-atom-types> <in-prevalence-file> <out-params-file> <out-pdbid-file> [--params=<params-file>] [--min-atom-types=<min-atom-types>] [--max-resolution=<max-resolution>] [--min-resolution=<min-resolution>] [--min-pdbids=<min-pdbids>] [--default-slope=<default-slope>]

Options:
    --params=<params-file>                  Overriding parameters file that includes radii, slopes, etc. [default: ]
    atom-type                               Generates atom types and other parameters from the chemical component descriptions.
    prevalence                              Generate prevalence of residues and atoms across a set of pdbids.
    parameters                              Generate parameters file.
    --residues=<comma-separated-residues>   Limit residues and atom-types to the given comma-separated residues. [default: ]
    --allow-errors                          Allow residues with errors.
    --min-atom-types=<min-atom-types>       Minimum number of a given atom-type required. [default: 50]
    --max-resolution=<max-resolution>       Maximum x-ray crystallographic resolution to allow. [default: 3.5]
    --min-resolution=<min-resolution>       Minimum x-ray crystallographic resolution to allow. [default: 0]
    --min-pdbids=<min-pdbids>               Minimum number of pdbids that have the minimum number of a given atom-type. [default: 1000]
    --default-slope=<default-slope>         Default b-factor correction slope to start with. [default: -0.5]

Unique atom-types are generated based on a "chemical coloring" approach that uses bond-type, element, estimated electrons, and aromatic information.

The wwPDB chemical components file is located at:
ftp://ftp.wwpdb.org/pub/pdb/data/monomers/components.cif
http://remediation.wwpdb.org/data/ccd

List of all wwPDB entry IDs along with entry type.
ftp://ftp.wwpdb.org/pub/pdb/derived_data/pdb_entry_type.txt
# Select only the protein, diffraction entries
% grep "diffraction" pdb_entry_type.txt | grep "prot" > pdb_prot_diffraction_only.txt

Amino Acid Residues: ALA,ARG,ASN,ASP,CYS,GLN,GLU,GLY,HIS,ILE,LEU,LYS,MET,PHE,PRO,SER,THR,TRP,TYR,VAL
Nucleic Acid Residues: A,C,G,I,U,DA,DC,DG,DI,DT,DU
Water: HOH
"""

import CifFile
import docopt
import urllib.request
import os
import collections
import json
import numpy as np

from . import __version__
from . import densityAnalysis

componentsFilename = "components.cif"

elementElectrons = {'H': 1, 'HE': 2, 'LI': 3, 'BE': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'NE': 10, 'NA': 11, 'MG': 12, 'AL': 13, 'SI': 14, 'P': 15, 'S': 16, 'CL': 17, 'AR': 18, 'K': 19, 'CA': 20, 'SC': 21, 'TI': 22,
                    'V': 23, 'CR': 24, 'MN': 25, 'FE': 26, 'CO': 27, 'NI': 28, 'CU': 29, 'ZN': 30, 'GA': 31, 'GE': 32, 'AS': 33, 'SE': 34, 'BR': 35, 'RB': 37, 'SR': 38, 'Y': 39, 'ZR': 40, 'NB': 41, 'MO': 42, 'TC': 43,
                    'RU': 44, 'RH': 45, 'PD': 46, 'AG': 47, 'CD': 48, 'IN': 49, 'SN': 50, 'SB': 51, 'TE': 52, 'I': 53, 'CS': 55, 'BA': 56, 'LA': 57, 'CE': 58, 'PR': 59, 'ND': 60, 'PM': 61, 'SM': 62, 'EU': 63, 'GD': 64,
                    'TB': 65, 'DY': 66, 'HO': 67, 'ER': 68, 'TM': 69, 'YB': 70, 'LU': 71, 'HF': 72, 'TA': 73, 'W': 74, 'RE': 75, 'OS': 76, 'IR': 77, 'PT': 78, 'AU': 79, 'HG': 80, 'TL': 81, 'PB': 82, 'BI': 83, 'PO': 84,
                    'RA': 88, 'AC': 89, 'TH': 90, 'PA': 91, 'U': 92, 'NP': 93, 'PU': 94, 'AM': 95}

elementAtomicRadii = {'H': 0.25, 'HE': 1.2, 'LI': 1.45, 'BE': 1.05, 'B': 0.85, 'C': 0.7, 'N': 0.65, 'O': 0.6, 'F': 0.5, 'NE': 1.6, 'NA': 1.8, 'MG': 1.5, 'AL': 1.25, 'SI': 1.1, 'P': 1.0, 'S': 1.0, 'CL': 1.0, 'AR': 0.71,
                      'K': 2.2, 'CA': 1.8, 'SC': 1.6, 'TI': 1.4, 'V': 1.35, 'CR': 1.4, 'MN': 1.4, 'FE': 1.4, 'CO': 1.35, 'NI': 1.35, 'CU': 1.35, 'ZN': 1.35, 'GA': 1.3, 'GE': 1.25, 'AS': 1.15, 'SE': 1.15, 'BR': 1.15,
                      'RB': 2.35, 'SR': 2.0, 'Y': 1.8, 'ZR': 1.55, 'NB': 1.45, 'MO': 1.45, 'TC': 1.35, 'RU': 1.3, 'RH': 1.35, 'PD': 1.4, 'AG': 1.6, 'CD': 1.55, 'IN': 1.55, 'SN': 1.45, 'SB': 1.45, 'TE': 1.4, 'I': 1.4,
                      'CS': 2.6, 'BA': 2.15, 'LA': 1.95, 'CE': 1.85, 'PR': 1.85, 'ND': 1.85, 'PM': 1.85, 'SM': 1.85, 'EU': 1.85, 'GD': 1.8, 'TB': 1.75, 'DY': 1.75, 'HO': 1.75, 'ER': 1.75, 'TM': 1.75, 'YB': 1.75,
                      'LU': 1.75, 'HF': 1.55, 'TA': 1.45, 'W': 1.35, 'RE': 1.35, 'OS': 1.3, 'IR': 1.35, 'PT': 1.35, 'AU': 1.35, 'HG': 1.5, 'TL': 1.9, 'PB': 1.8, 'BI': 1.6, 'PO': 1.9, 'RA': 2.15, 'AC': 1.95, 'TH': 1.8,
                      'PA': 1.8, 'U': 1.75, 'NP': 1.75, 'PU': 1.75, 'AM': 1.75}

oxygenDoubleBondColor = "O.N.8.DOUB.N"
oxygenSingleBondColor = "O.N.9.SING.N"

def main():
    args = docopt.docopt(__doc__, version=__version__)
    if args["--help"]:
        print(__doc__)
        exit(0)

    if args["atom-type"]:
        if os.path.isfile("components_info.json"):
            with open("components_info.json", 'r') as jsonFile:
                componentsInfo = json.load(jsonFile)
        else:
            componentsInfo = processComponents()

        if not os.path.isfile("components_info.json"):
            with open("components_info.json", 'w') as outFile:
                print(json.dumps(componentsInfo, indent=2, sort_keys=True), file=outFile)

        componentsInfo["errors"] = set(componentsInfo["errors"])
        allowedResidueTypes = set(args["--residues"].split(",")) if args["--residues"] else set()
        initialParameters = { "residue_name_map_electrons" : {}, "full_atom_name_map_atom_type" : {}, "full_atom_name_map_electrons" : {}, "element_map_electrons" : elementElectrons, "leaving_atoms" : [] }

        for residue in componentsInfo["residues"].values():
            residue["estimated_electrons"] = 0
            for atom in residue["atoms"].values():
                atom["num_bound_hydrogens"] = sum(1 for atomName,bondType,aromatic,stereo in atom["bonds"] if atomName in residue["atoms"] and residue["atoms"][atomName]["element"] == "H" and residue["atoms"][atomName]["leaving"] == atom["leaving"])
                if atom["element"] in elementElectrons:
                    try:
                        atom["charge"] = float(atom["charge"])
                    except:
                        atom["charge"] = 0
                    atom["estimated_electrons"] = elementElectrons[atom["element"]] + atom["num_bound_hydrogens"] - atom["charge"]
                else:
                    atom["estimated_electrons"] = 0
                    componentsInfo["errors"].add(residue["name"])

            for atom in residue["atoms"].values():
                atom["element_color"] = atom["element"] + "." + atom["aromatic"] + "." + str(int(atom["estimated_electrons"]))
            for atom in residue["atoms"].values():
                atom["element_bond_colors"] = [ residue["atoms"][atomName]["element_color"] + "." +  bondType + "." + aromatic for atomName,bondType,aromatic,stereo in atom["bonds"]
                                                if atomName in residue["atoms"] and (atom["leaving"] == "Y" or atom["leaving"] == residue["atoms"][atomName]["leaving"]) ]
                atom["full_element_color"] = atom["element_color"] + "#" + "_".join(sorted(atom["element_bond_colors"]))

            # identify and adjust oxygen atoms in resonance.
            for testAtom in residue["atoms"].values():
                if oxygenDoubleBondColor in testAtom["element_bond_colors"] and oxygenSingleBondColor in testAtom["element_bond_colors"]:
                    oxygenAtomTuples = [(residue["atoms"][atomName],residue["atoms"][atomName]["element_color"] + "." +  bondType + "." + aromatic) for atomName,bondType,aromatic,stereo in testAtom["bonds"]
                                        if atomName in residue["atoms"] and residue["atoms"][atomName]["element"] == "O"  and testAtom["leaving"] == "Y" or testAtom["leaving"] == residue["atoms"][atomName]["leaving"] ]
                    resonanceOxygenAtoms = [atom for atom,bondColor in oxygenAtomTuples if bondColor == oxygenDoubleBondColor or bondColor == oxygenSingleBondColor ]
                    if len(set(atom["estimated_electrons"] for atom in resonanceOxygenAtoms)) > 1:
                        averageElectrons = np.mean([atom["estimated_electrons"] for atom in resonanceOxygenAtoms])
                        longestFullElementColor = ""
                        for atom in resonanceOxygenAtoms:
                            atom["estimated_electrons"] = averageElectrons
                            atom["element_color"] = atom["element"] + "." + atom["aromatic"] + "." + str(float(atom["estimated_electrons"]))[:5]
                            atom["element_bond_colors"] = [ residue["atoms"][atomName]["element_color"] + "." +  "RESON" + "." + aromatic for atomName,bondType,aromatic,stereo in atom["bonds"]
                                                            if atomName in residue["atoms"]  and (atom["leaving"] == "Y" or atom["leaving"] == residue["atoms"][atomName]["leaving"]) ]
                            atom["full_element_color"] = atom["element_color"] + "#" + "_".join(sorted(atom["element_bond_colors"]))
                            if len(atom["full_element_color"]) > len(longestFullElementColor):
                                longestFullElementColor = atom["full_element_color"]
                        for atom in resonanceOxygenAtoms:
                            atom["full_element_color"] = longestFullElementColor

            if args["--allow-errors"] or residue["name"] not in componentsInfo["errors"]:
                for atom in residue["atoms"].values():
                    residue["estimated_electrons"] += atom["estimated_electrons"] if atom["leaving"] != "Y" else 0
                residue["estimated_electrons"] = float(np.round(residue["estimated_electrons"]))

                initialParameters["residue_name_map_electrons"][residue["name"]] = residue["estimated_electrons"]
                if not allowedResidueTypes or residue["name"] in allowedResidueTypes: # only add certain residues
                    for atom in residue["atoms"].values():
                        if atom["element"] != "H":
                            fullAtomName = residue["name"] + "_" + atom["name"]
                            initialParameters["full_atom_name_map_atom_type"][fullAtomName] = atom["full_element_color"]
                            initialParameters["full_atom_name_map_electrons"][fullAtomName] = atom["estimated_electrons"]
                            if atom["leaving"] == "Y":
                                initialParameters["leaving_atoms"].append(fullAtomName)

        # Print report
        print("Unique Residue Types:",len(set(name.split("_")[0] for name in initialParameters["full_atom_name_map_atom_type"].keys())))
        print("Unique Full Atom Names:",len(set(initialParameters["full_atom_name_map_atom_type"].keys())))
        print("Unique Atom Types:",len(set(initialParameters["full_atom_name_map_atom_type"].values())))
        # testElectrons = collections.defaultdict(set)
        # for name,value in initialParameters["full_atom_name_map_atom_type"].items():
        #     testElectrons[value].add(initialParameters["full_atom_name_map_electrons"][name])
        # for atomType,values in testElectrons.items():
        #     if len(values) > 1:
        #         print(atomType,values)

        with open(args["<out-jsonfile>"], 'w') if args["<out-jsonfile>"] != "-" else sys.stdout as outFile:
            print(json.dumps(initialParameters, indent=2, sort_keys=True), file=outFile)
    elif args["prevalence"]:
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
    elif args["parameters"]:
        if args["--params"]:
            try:
                with open(args["--params"], 'r') as paramsFile:
                    overrideParams = json.load(paramsFile)
            except:
                sys.exit(str("Error: params file \"") + args["--params"] + "\" does not exist or is not parsable.")

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
        residueName = residueName.upper()
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



