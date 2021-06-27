"""
The generate parameters mode command line interface
  Generates initial parameters needed for pdb_eda calculations.

Usage:
    pdb_eda generate -h | --help
    pdb_eda generate atom-type <out-jsonfile> [--residues=<comma-separated-residues>] [--allow-errors] [--default-slope=<default-slope>] [--F000]
    pdb_eda generate prevalence <pdbid-file> <out-jsonfile> [--testing]
    pdb_eda generate parameters <in-atom-types> <in-prevalence-file> <out-params-file> <out-pdbid-file> [--params=<params-file>] [--min-atom-types=<min-atom-types>] [--min-atoms=<min-atoms>] [--max-atoms=<max-atoms>] [--max-resolution=<max-resolution>] [--min-resolution=<min-resolution>]

Options:
    atom-type                               Generates atom types and other parameters from the chemical component descriptions.
    --residues=<comma-separated-residues>   Limit residues and atom-types to the given comma-separated residues. [default: ]
    --allow-errors                          Allow residues with errors.
    --default-slope=<default-slope>         Default b-factor correction slope to start with. [default: -0.5]
    --F000                                  Generate parameters for F000 calculations.
    prevalence                              Generate prevalence report of residues and atoms across a set of pdbids.
    --testing                               Run only a single process for testing purposes.
    parameters                              Generate parameters file.
    --params=<params-file>                  Overriding parameters file that includes radii, slopes, etc. Useful for including currently optimized atom types. [default: ]
    --min-atom-types=<min-atom-types>       Minimum number of a given atom-type required. [default: 5]
    --min-atoms=<min-atoms>                 Minimum number of parameter atoms required. [default: 500]
    --max-atoms=<max-atoms>                 Maximum number of (all) parameter atoms allowed. Useful to selecting entries that are reasonably fast to analyze [default: 5000]
    --max-resolution=<max-resolution>       Maximum x-ray crystallographic resolution to allow. [default: 3.5]
    --min-resolution=<min-resolution>       Minimum x-ray crystallographic resolution to allow. [default: 0]


Unique atom-types are generated based on a "chemical coloring" approach that uses bond-type, element, estimated electrons, and aromatic information.

The wwPDB chemical components file is located at:
ftp://ftp.wwpdb.org/pub/pdb/data/monomers/components.cif
http://remediation.wwpdb.org/data/ccd

List of all wwPDB entry IDs along with entry type.
ftp://ftp.wwpdb.org/pub/pdb/derived_data/pdb_entry_type.txt
# Select only the protein, diffraction entries
% grep "diffraction" pdb_entry_type.txt | grep "prot" > pdb_prot_diffraction_only.txt

# Faster to download all gzipped compressed PDB entries directly from wwPDB before analyzing
% cd pdb_data
% rsync -rLpt -v -z --port=33444 rsync.rcsb.org::ftp_data/structures/all/pdb/ ./

Amino Acid Residues: ALA,ARG,ASN,ASP,CYS,GLN,GLU,GLY,HIS,ILE,LEU,LYS,MET,PHE,PRO,SER,THR,TRP,TYR,VAL
Nucleic Acid Residues: A,C,G,I,U,DA,DC,DG,DI,DT,DU

Typical order of use:
1) atom-type: generate new set of atom types to optimize.
2) prevalence: analyze a set of PDB entries for atom-type prevalence.
3) parameters: generate new parameters file and list of PDB IDs for optimization.
"""

import CifFile
import docopt
import urllib.request
import os
import collections
import json
import numpy as np
import multiprocessing

from . import __version__
from . import densityAnalysis
from . import fileUtils

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

oxygenDoubleBondColor = "O.N.8.DOUB"
oxygenSingleBondColor = "O.N.9.SING"

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
        args["--default-slope"] = float(args["--default-slope"])

        if args["--F000"]:
            initialParams = {"full_atom_name_map_electrons": {}, "element_map_electrons": elementElectrons }
        else:
            initialParams = { "full_atom_name_map_atom_type" : {}, "full_atom_name_map_electrons" : {}, "leaving_atoms" : [], "radii" : {}, "slopes" : {}, "bonded_atoms" : {} }

        fullAtomName2AtomType = {}

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
                atom["element_bond_colors"] = [ residue["atoms"][atomName]["element_color"] + "." +  bondTyping(bondType,aromatic) for atomName,bondType,aromatic,stereo in atom["bonds"]
                                                if atomName in residue["atoms"] and (atom["leaving"] == "Y" or atom["leaving"] == residue["atoms"][atomName]["leaving"]) ]
                atom["full_element_color"] = atom["element_color"] + "#" + "_".join(sorted(atom["element_bond_colors"]))

            # identify and adjust oxygen atoms in resonance.
            for testAtom in residue["atoms"].values():
                if oxygenDoubleBondColor in testAtom["element_bond_colors"] and oxygenSingleBondColor in testAtom["element_bond_colors"]:
                    oxygenAtomTuples = [(residue["atoms"][atomName],residue["atoms"][atomName]["element_color"] + "." + bondTyping(bondType,aromatic)) for atomName,bondType,aromatic,stereo in testAtom["bonds"]
                                        if atomName in residue["atoms"] and residue["atoms"][atomName]["element"] == "O"  and testAtom["leaving"] == "Y" or testAtom["leaving"] == residue["atoms"][atomName]["leaving"] ]
                    resonanceOxygenAtoms = [atom for atom,bondColor in oxygenAtomTuples if bondColor == oxygenDoubleBondColor or bondColor == oxygenSingleBondColor ]
                    if len(set(atom["estimated_electrons"] for atom in resonanceOxygenAtoms)) > 1:
                        averageElectrons = np.mean([atom["estimated_electrons"] for atom in resonanceOxygenAtoms])
                        longestFullElementColor = ""
                        for atom in resonanceOxygenAtoms:
                            atom["estimated_electrons"] = averageElectrons
                            atom["element_color"] = atom["element"] + "." + atom["aromatic"] + "." + str(float(atom["estimated_electrons"]))[:5]
                            atom["element_bond_colors"] = [ residue["atoms"][atomName]["element_color"] + "." +  "RESON" for atomName,bondType,aromatic,stereo in atom["bonds"]
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

                if not allowedResidueTypes or residue["name"] in allowedResidueTypes: # only add certain residues
                    for atom in residue["atoms"].values():
                        if atom["element"] != "H":
                            fullAtomName = residue["name"].strip() + "_" + atom["name"]
                            initialParams["full_atom_name_map_electrons"][fullAtomName] = atom["estimated_electrons"]
                            fullAtomName2AtomType[fullAtomName] = atom["full_element_color"]
                            if not args["--F000"]:
                                initialParams["full_atom_name_map_atom_type"][fullAtomName] = atom["full_element_color"]
                                initialParams["radii"][atom["full_element_color"]] = elementAtomicRadii[atom["element"]]
                                initialParams["slopes"][atom["full_element_color"]] = args["--default-slope"]
                                initialParams["bonded_atoms"][fullAtomName] = [residue["name"].strip() + "_" + atomName for atomName, bondType, aromatic, stereo in atom["bonds"] if residue["atoms"][atomName]["element"] != "H"]
                                if atom["leaving"] == "Y":
                                    initialParams["leaving_atoms"].append(fullAtomName)

        # Print report
        print("Unique Residue Types:",len(set(name.split("_")[0] for name in fullAtomName2AtomType.keys())))
        print("Unique Full Atom Names:",len(set(fullAtomName2AtomType.keys())))
        print("Unique Atom Types:",len(set(fullAtomName2AtomType.values())))
        # testElectrons = collections.defaultdict(set)
        # for name,value in initialParams["full_atom_name_map_atom_type"].items():
        #     testElectrons[value].add(initialParams["full_atom_name_map_electrons"][name])
        # for atomType,values in testElectrons.items():
        #     if len(values) > 1:
        #         print(atomType,values)

        with open(args["<out-jsonfile>"], 'w') if args["<out-jsonfile>"] != "-" else sys.stdout as outFile:
            print(json.dumps(initialParams, indent=2, sort_keys=True), file=outFile)
    elif args["prevalence"]:
        try:
            pdbids = []
            with open(args["<pdbid-file>"], "r") as textFile:
                for pdbid in textFile:
                    pdbids.append(pdbid[0:4])
        except:
            RuntimeError(str("Error: PDB IDs file \"") + args["<pdbid-file>"] + "\" does not exist or is not parsable.")


        if args["--testing"]:
            results = [ processFunction(pdbid) for pdbid in pdbids ]
        else:
            with multiprocessing.Pool() as pool:
                results = pool.map(processFunction, pdbids)

        pdbidInfo = {}
        for resultFilename in results:
            if resultFilename:
                try:
                    with open(resultFilename, 'r') as jsonFile:
                        result = json.load(jsonFile)
                        pdbidInfo[result["pdbid"]] = result
                    os.remove(resultFilename)
                except:
                    pass

        totalResidueAtoms = collections.defaultdict(int)
        totalElements = collections.defaultdict(int)
        totalResidues = collections.defaultdict(int)
        for pdbid,info in pdbidInfo.items():
            for name,count in info["full_atom_name_counts"].items():
                totalResidueAtoms[name] += count
            for name,count in info["element_counts"].items():
                totalElements[name] += count
            for name,count in info["residue_counts"].items():
                totalResidues[name] += count

        with open(args["<out-jsonfile>"], 'w') if args["<out-jsonfile>"] != "-" else sys.stdout as outFile:
            jsonOutput = { "pdbid_info" : pdbidInfo, "full_atom_name_counts" : totalResidueAtoms, "residue_counts" : totalResidues, "element_counts" : totalElements }
            print(json.dumps(jsonOutput, indent=2, sort_keys=True), file=outFile)
    elif args["parameters"]:
        args["--max-resolution"] = float(args["--max-resolution"])
        args["--min-resolution"] = float(args["--min-resolution"])
        args["--min-atom-types"] = int(args["--min-atom-types"])
        args["--min-atoms"] = int(args["--min-atoms"])
        args["--max-atoms"] = int(args["--max-atoms"])

        if args["--params"]:
            try:
                with open(args["--params"], 'r') as paramsFile:
                    overrideParams = json.load(paramsFile)
            except:
                RuntimeError(str("Error: params file \"") + args["--params"] + "\" does not exist or is not parsable.")
        else:
            overrideParams = None

        try:
            with open(args["<in-atom-types>"], 'r') as paramsFile:
                initialParams = json.load(paramsFile)
        except:
            RuntimeError(str("Error: params file \"") + args["<in-atom-types>"] + "\" does not exist or is not parsable.")

        try:
            with open(args["<in-prevalence-file>"], 'r') as jsonFile:
                prevalenceInfo = json.load(jsonFile)
        except:
            RuntimeError(str("Error: prevalence file \"") + args["<in-prevalence-file>"] + "\" does not exist or is not parsable.")

        currentPDBinfo = { pdbid:info for pdbid,info in prevalenceInfo["pdbid_info"].items()
                           if info["properties"]["resolution"] >= args["--min-resolution"] and info["properties"]["resolution"] <= args["--max-resolution"] }

        testingFullAtomNames = [fullAtomName for fullAtomName in initialParams["full_atom_name_map_atom_type"].keys()
                                if fullAtomName not in initialParams["leaving_atoms"] and (not overrideParams or fullAtomName not in overrideParams["full_atom_name_map_atom_type"])]
        testingAtomTypes = set(initialParams["full_atom_name_map_atom_type"][fullAtomName] for fullAtomName in testingFullAtomNames)

        allFullAtomNames = list(testingFullAtomNames)
        if overrideParams:
            allFullAtomNames.extend([fullAtomName for fullAtomName in overrideParams["full_atom_name_map_atom_type"].keys() if fullAtomName not in overrideParams["leaving_atoms"]])

        pdbids = []
        for pdbid, info in currentPDBinfo.items():
            atomTypeSum = { atomType:0 for atomType in testingAtomTypes }
            for fullAtomName in testingFullAtomNames:
                atomTypeSum[initialParams["full_atom_name_map_atom_type"][fullAtomName]] += info["full_atom_name_counts"][fullAtomName] if fullAtomName in info["full_atom_name_counts"] else 0

            analyzableAtoms = sum(atomTypeSum.values())
            totalAtoms = sum(info["full_atom_name_counts"][fullAtomName] for fullAtomName in allFullAtomNames if fullAtomName in info["full_atom_name_counts"])
            if all(count >= args["--min-atom-types"] for count in atomTypeSum.values()) and analyzableAtoms >= args["--min-atoms"] and totalAtoms <= args["--max-atoms"]:
                pdbids.append(pdbid)

        with open(args['<out-pdbid-file>'], "w") if args["<out-pdbid-file>"] != "-" else sys.stdout as txtFile:
            print("\n".join(pdbids), file=txtFile)

        if overrideParams:
            initialParams["full_atom_name_map_atom_type"].update(overrideParams["full_atom_name_map_atom_type"])
            initialParams["full_atom_name_map_electrons"].update(overrideParams["full_atom_name_map_electrons"])
            initialParams["radii"].update(overrideParams["radii"])
            initialParams["slopes"].update(overrideParams["slopes"])
            leavingAtoms = set(initialParams["leaving_atoms"])
            leavingAtoms.update(overrideParams["leaving_atoms"])
            initialParams["leaving_atoms"] = list(leavingAtoms)
            initialParams["optimize"] = [atomType for atomType in initialParams["radii"].keys() if atomType not in overrideParams["radii"]]

        with open(args["<out-params-file>"], 'w') if args["<out-params-file>"] != "-" else sys.stdout as outFile:
            print(json.dumps(initialParams, indent=2, sort_keys=True), file=outFile)


def processComponents():
    """Retrieves and extracts residues, atoms, and bonds from the wwPDB chemical components Cif file.

    :return: component information
    :rtype: :py:class:`dict`
    """
    if not os.path.isfile(componentsFilename):
        urllib.request.urlretrieve("ftp://ftp.wwpdb.org/pub/pdb/data/monomers/components.cif", componentsFilename)

    print("Parsing the components.cif file.  This takes a few minutes.")
    components = CifFile.ReadCif(componentsFilename)

    print("Processing the components.")
    residues = {}
    errors = set()
    for residueName in components.keys():
        residue = components[residueName]
        residueName = residueName.upper().strip()
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


def bondTyping(bondType, aromatic):
    """Returns type of bond based on given bondType and whether it was marked aromatic.

    :param bondType:
    :type bondType: :py:class:`str`
    :param aromatic:
    :type aromatic: :py:class:`str`

    :return: bondType
    :rtype: :py:class:`str`
    """
    return bondType if aromatic == "N" else "AROM"

def processFunction(pdbid):
    """Process function for multiprocessing to analyze a single pdb entry.

    :param pdbid: pdbid for entry to download and analyze.
    :type pdbid: :py:class:`str`

    :return: resultFilename
    :rtype: :py:class:`str`
    """
    if not densityAnalysis.testCCP4URL(pdbid):
        return 0

    analyzer = densityAnalysis.fromPDBid(pdbid, ccp4density=False, ccp4diff=False)
    if not analyzer:
        return 0

    info = {}
    info["pdbid"] = pdbid
    info["properties"] = {property: value for (property, value) in analyzer.biopdbObj.header.items()}
    info["properties"]["resolution"] = float(analyzer.pdbObj.header.resolution)
    info["properties"]['space_group'] = analyzer.pdbObj.header.spaceGroup
    info["full_atom_name_counts"] = collections.Counter(densityAnalysis.residueAtomName(atom) for residue in analyzer.biopdbObj.get_residues() for atom in residue.get_atoms())
    info["element_counts"] = collections.Counter(atom.element for residue in analyzer.biopdbObj.get_residues() for atom in residue.get_atoms())
    info["residue_counts"] = collections.Counter(residue.resname for residue in analyzer.biopdbObj.get_residues())

    resultFilename = fileUtils.createTempJSONFile(info, "tempResults_")
    return resultFilename
