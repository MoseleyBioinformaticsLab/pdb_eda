"""
pdb_eda multiple structure analysis mode command-line interface
  Analyzes multiple pdb entries.

Usage:
    pdb_eda multiple -h | --help
    pdb_eda multiple <pdbid-file> <out-result-file> [--params=<params-file>] [--out-format=<format>] [--testing] [--time-out=<seconds>] [--silent]
    pdb_eda multiple <in-result-file> <out-pdbid-file> --filter [--out-format=<format>] [--max-resolution=<max-resolution>] [--min-resolution=<min-resolution>] [--min-atoms=<min-atoms>] [--min-residues=<min-residues>] [--min-elements=<min-elements>]
    pdb_eda multiple <pdbid-file> --reload
    pdb_eda multiple <pdbid-file> <out-dir> --single-mode=<quoted-single-mode-options> [--time-out=<seconds>] [--testing] [--silent]
    pdb_eda multiple <pdbid-file> <out-dir> --contacts-mode=<quoted-contacts-mode-options> [--time-out=<seconds>] [--testing] [--silent]

Options:
    -h, --help                                      Show this screen.
    <out-result-file>                               Output filename. "-" will write to standard output.
    <in-result-file>                                Input filename. "-" will read from standard input.  Format must be JSON.
    <pdbid-file>                                    Input filename that contains the pdb ids. "-" will read from standard input.
    <out-pdbid-file>                                Output filename that contains the pdb ids. "-" will write to standard output.
    --params=<params-file>                          Overriding parameters file that includes radii, slopes, etc. [default: ]
    --out-format=<format>                           Output file format, available formats: csv, json [default: json].
    --time-out=<seconds>                            Set a maximum time to try to analyze any single pdb entry. [default: 0]
    --testing                                       Run only a single process for testing purposes.
    --silent                                        Do not print error message to stderr.
    --single-mode=<quoted-single-mode-options>      Run single structure analysis mode on a set of PDB entries.
    --contacts-mode=<quoted-contacts-mode-options>  Run (crystal) contacts analysis mode on a set of PDB entries.
    --filter                                        Filter pdbids based on criteria.
    --max-resolution=<max-resolution>               Maximum x-ray crystallographic resolution to allow. [default: 3.5]
    --min-resolution=<min-resolution>               Minimum x-ray crystallographic resolution to allow. [default: 0]
    --min-atoms=<min-atoms>                         Minimum number of atoms analyzed (numAtomsAnalyzed) required. [default: 300]
    --min-residues=<min-residues>                   Minimum number of residues required of the given residue types separated by commas. [default: 0]
    --min-elements=<min-elements>                   Minimum number of atoms required of the given elements separated by commas. [default: 0]
    --reload                                        Test that each pdbid loads and try to download again any pdbid that fail.

List of all wwPDB entry IDs along with entry type.
ftp://ftp.wwpdb.org/pub/pdb/derived_data/pdb_entry_type.txt
"""

import os
import gc
import time
import sys

import json
import csv
import multiprocessing
import signal
import collections
import docopt

from . import densityAnalysis
from . import __version__
from . import singleStructure
from . import crystalContacts
from . import fileUtils

globalArgs = {}

def main():
    global globalArgs
    globalArgs = docopt.docopt(__doc__, version=__version__)
    if globalArgs["--help"]:
        print(__doc__)
        exit(0)
    globalArgs["--time-out"] = int(globalArgs["--time-out"])

    if globalArgs["--params"]:
        try:
            with open(globalArgs["--params"], 'r') as paramsFile:
                params = json.load(paramsFile)
            densityAnalysis.setGlobals(params)
        except:
            raise RuntimeError(str("Error: params file \"") + globalArgs["--params"] + "\" does not exist or is not parsable.")

    if globalArgs["--filter"]:
        globalArgs["--max-resolution"] = float(globalArgs["--max-resolution"])
        globalArgs["--min-resolution"] = float(globalArgs["--min-resolution"])
        globalArgs["--min-atoms"] = int(globalArgs["--min-atoms"])

        if "," in globalArgs["--min-residues"]:
            allowedResidueTypes = globalArgs["--min-residues"].split(",")
            globalArgs["--min-residues"] = float(allowedResidueTypes.pop(0))
            allowedResidueTypes = set(allowedResidueTypes)
        else:
            globalArgs["--min-residues"] = float(globalArgs["--min-residues"])
            allowedResidueTypes = set()

        if "," in globalArgs["--min-elements"]:
            allowedElements = globalArgs["--min-elements"].split(",")
            globalArgs["--min-elements"] = float(allowedElements.pop(0))
            allowedElements = set(allowedElements)
        else:
            globalArgs["--min-elements"] = float(globalArgs["--min-elements"])
            allowedElements = set()

        try:
            with open(globalArgs["<in-result-file>"], "r") if globalArgs["<in-result-file>"] != "-" else sys.stdin as jsonFile:
                pdbResults = json.load(jsonFile)
        except:
            raise RuntimeError(str("Error: PDB results file \"") + globalArgs["<in-result-file>"] + "\" does not exist or is not parsable.")

        pdbids = [pdbid for pdbid in pdbResults if
                  pdbResults[pdbid]["stats"]["num_atoms_analyzed"] >= globalArgs["--min-atoms"] and
                  float(pdbResults[pdbid]["stats"]["resolution"]) >= globalArgs["--min-resolution"] and
                  float(pdbResults[pdbid]["stats"]["resolution"]) <= globalArgs["--max-resolution"] and
                  sum(count for (residueType, count) in pdbResults[pdbid]["properties"]["residue_counts"].items() if not allowedResidueTypes or residueType in allowedResidueTypes) >= globalArgs["--min-residues"] and
                  sum(count for (element, count) in pdbResults[pdbid]["properties"]["element_counts"].items() if not allowedElements or element in allowedElements) >= globalArgs["--min-elements"]
                  ]

        if globalArgs["--out-format"] == "json":
            with open(globalArgs['<out-pdbid-file>'], "w") if globalArgs["<out-pdbid-file>"] != "-" else sys.stdout as jsonFile:
                print(json.dumps(pdbids, indent=2, sort_keys=True), file=jsonFile)
        else:
            with open(globalArgs['<out-pdbid-file>'], "w") if globalArgs["<out-pdbid-file>"] != "-" else sys.stdout as txtFile:
                print("\n".join(pdbids), file=txtFile)
    elif globalArgs["--reload"]:
        try:
            pdbids = []
            with open(globalArgs["<pdbid-file>"], "r") if globalArgs["<pdbid-file>"] != "-" else sys.stdin as textFile:
                for pdbid in textFile:
                    pdbids.append(pdbid[0:4])
        except:
            raise RuntimeError(str("Error: PDB IDs file \"") + globalArgs["<pdbid-file>"] + "\" does not exist or is not parsable.")

        badPDBids = [pdbid for pdbid in pdbids if not testPDBIDLoad(pdbid)]
        for pdbid in badPDBids:
            densityAnalysis.cleanPDBid(pdbid)
        badPDBids = [pdbid for pdbid in badPDBids if not testPDBIDLoad(pdbid)]
        for pdbid in badPDBids:
            densityAnalysis.cleanPDBid(pdbid)
        if badPDBids:
            print("Bad PDBids:",",".join(badPDBids))
    else:
        try:
            pdbids = []
            with open(globalArgs["<pdbid-file>"], "r") if globalArgs["<pdbid-file>"] != "-" else sys.stdin as textFile:
                for pdbid in textFile:
                    pdbids.append(pdbid[0:4])
        except:
            raise RuntimeError(str("Error: PDB IDs file \"") + globalArgs["<pdbid-file>"] + "\" does not exist or is not parsable.")

        if globalArgs["--single-mode"]:
            processFunction = singleModeFunction
            if not os.path.isdir(globalArgs["<out-dir>"]):
                if not os.path.isfile(globalArgs["<out-dir>"]):
                    os.mkdir(globalArgs["<out-dir>"])
                else:
                    raise RuntimeError(str("Error: Output directory \"") + globalArgs["<out-dir>"] + "\" is a file.")
        elif globalArgs["--contacts-mode"]:
            processFunction = contactsModeFunction
            if not os.path.isdir(globalArgs["<out-dir>"]):
                if not os.path.isfile(globalArgs["<out-dir>"]):
                    os.mkdir(globalArgs["<out-dir>"])
                else:
                    raise RuntimeError(str("Error: Output directory \"") + globalArgs["<out-dir>"] + "\" is a file.")
        else:
            processFunction = multipleModeFunction

        if globalArgs["--testing"]:
            results = [ processFunction(pdbid) for pdbid in pdbids ]
        else:
            with multiprocessing.Pool() as pool:
                results = pool.map(processFunction, pdbids, chunksize=1)

        if not globalArgs["--single-mode"] and not globalArgs["--contacts-mode"]: # skip generating results if running in another mode or submode.
            fullResults = {}
            for resultFilename in results:
                if resultFilename:
                    try:
                        with open(resultFilename, 'r') as jsonFile:
                            result = json.load(jsonFile)
                            fullResults[result["pdbid"]] = result
                        os.remove(resultFilename)
                    except:
                        pass

            if globalArgs["--out-format"] == 'csv' or globalArgs["--out-format"] == 'txt':
                statsHeaders = ['density_electron_ratio', 'voxel_volume', 'f000', 'num_voxels_aggregated', 'total_aggregated_electrons', 'density_mean', 'diff_density_mean', 'resolution', 'space_group',
                                'num_atoms_analyzed', 'num_residue_clouds_analyzed', 'num_domain_clouds_analyzed', "atom_overlap_completeness"]
                with open(globalArgs['<out-result-file>'], "w", newline='') if globalArgs["<out-result-file>"] != "-" else sys.stdout as csvFile:
                    writer = csv.writer(csvFile)
                    writer.writerow(['pdbid'] + statsHeaders + sorted(densityAnalysis.paramsGlobal["radii"]))
                    for result in fullResults.values():
                        stats = [result["stats"][header] for header in statsHeaders]
                        diffs = [result["diffs"][atomType] for atomType in sorted(densityAnalysis.paramsGlobal["radii"])]
                        writer.writerow([result['pdbid']] + stats + diffs)
            else:
                with open(globalArgs['<out-result-file>'], "w") if globalArgs["<out-result-file>"] != "-" else sys.stdout as jsonFile:
                    print(json.dumps(fullResults, indent=2, sort_keys=True), file=jsonFile)


def singleModeFunction(pdbid):
    """Execute pdb_eda single mode on a given pdbid.

    :param pdbid:
    :type pdbid: :py:class:`str`

    :return: 0
    """
    singleModeCommandLine = "pdb_eda single " + pdbid + " " + globalArgs["<out-dir>"] + "/" + pdbid + ".result " + globalArgs["--single-mode"]
    sys.argv = singleModeCommandLine.split()
    if globalArgs["--time-out"]:
        try:
            with timeout(seconds=globalArgs["--time-out"]):
                singleStructure.main()
        except Exception as exception:
            if not globalArgs["--silent"]:
                print(pdbid, exception, file=sys.stderr)
    elif globalArgs["--testing"]:
        singleStructure.main()
    else:
        try:
            singleStructure.main()
        except Exception as exception:
            if not globalArgs["--silent"]:
                print(pdbid, exception, file=sys.stderr)

    gc.collect()
    return 0


def contactsModeFunction(pdbid):
    """Execute pdb_eda contacts mode on a given pdbid.

    :param pdbid:
    :type pdbid: :py:class:`str`

    :return: 0
    """
    contactsModeCommandLine = "pdb_eda contacts " + pdbid + " " + globalArgs["<out-dir>"] + "/" + pdbid + ".result " + globalArgs["--contacts-mode"]
    sys.argv = contactsModeCommandLine.split()
    if globalArgs["--time-out"]:
        try:
            with timeout(seconds=globalArgs["--time-out"]):
                crystalContacts.main()
        except Exception as exception:
            if not globalArgs["--silent"]:
                print(pdbid, exception, file=sys.stderr)
    elif globalArgs["--testing"]:
        crystalContacts.main()
    else:
        try:
            crystalContacts.main()
        except Exception as exception:
            if not globalArgs["--silent"]:
                print(pdbid, exception, file=sys.stderr)

    gc.collect()
    return 0


def multipleModeFunction(pdbid):
    """Analyze a single pdbid.

    :param pdbid:
    :type pdbid: py:class:`str`

    :return: resultFilename or 0
    :rtype: :py:class:`str`, :py:class:`int`
    """
    if globalArgs["--time-out"]:
        try:
            with timeout(seconds=globalArgs["--time-out"]):
                return analyzePDBID(pdbid)
        except Exception as exception:
            if not globalArgs["--silent"]:
                print(pdbid, exception, file=sys.stderr)
            return 0
    else:
        return analyzePDBID(pdbid)

def testPDBIDLoad(pdbid):
    """Returns whether the pdbid downloads and loads correctly.

    :param pdbid:
    :type pdbid: :py:class:`str`

    :return: bool
    :rtype: :py:class:`bool`
    """
    analyzer = densityAnalysis.fromPDBid(pdbid)
    return False if not analyzer else True

def analyzePDBID(pdbid):
    """Process function to analyze a single pdb entry.

    :param pdbid: pdbid for entry to download and analyze.
    :type pdbid: :py:class:`str`

    :return: resultFilename
    :rtype: :py:class:`str`
    """
    startTime = time.process_time()

    analyzer = densityAnalysis.fromPDBid(pdbid)
    if not analyzer or not analyzer.densityElectronRatio:
        return 0

    diffs = {atomType:((analyzer.medians['corrected_density_electron_ratio'][atomType] - analyzer.densityElectronRatio) / analyzer.densityElectronRatio)
              if atomType in analyzer.medians['corrected_density_electron_ratio'] else 0 for atomType in sorted(densityAnalysis.paramsGlobal["radii"])}

    atomOverlapCompleteness = sum(analyzer.atomTypeOverlapCompleteness.values())
    atomOverlapInCompleteness = sum(analyzer.atomTypeOverlapIncompleteness.values())
    if atomOverlapCompleteness > 0 or atomOverlapInCompleteness > 0:
        atomOverlapCompleteness = atomOverlapCompleteness / (atomOverlapCompleteness + atomOverlapInCompleteness)

    stats = {'density_electron_ratio' : analyzer.densityElectronRatio, 'voxel_volume' : analyzer.densityObj.header.unitVolume, 'f000' : analyzer.F000, 'num_voxels_aggregated' : analyzer.numVoxelsAggregated,
        'total_aggregated_electrons' : analyzer.totalAggregatedElectrons, 'density_mean' : analyzer.densityObj.header.densityMean, 'diff_density_mean' : analyzer.diffDensityObj.header.densityMean,
        'resolution' : analyzer.pdbObj.header.resolution, 'space_group' : analyzer.pdbObj.header.spaceGroup, 'num_atoms_analyzed' : len(analyzer.atomCloudDescriptions),
        'num_residue_clouds_analyzed' : len(analyzer.residueCloudDescriptions), 'num_domain_clouds_analyzed' : len(analyzer.domainCloudDescriptions), 'atom_overlap_completeness' : atomOverlapCompleteness}

    properties = { property : value for (property,value) in analyzer.biopdbObj.header.items() }
    properties['residue_counts'] = dict(collections.Counter(residue.resname for residue in analyzer.biopdbObj.get_residues()))
    properties['element_counts'] = dict(collections.Counter(atom.element for atom in analyzer.biopdbObj.get_atoms()))

    elapsedTime = time.process_time() - startTime
    resultFilename = fileUtils.createTempJSONFile({ "pdbid" : analyzer.pdbid, "diffs" : diffs, "stats" : stats, "execution_time" : elapsedTime, "properties" : properties }, "tempResults_")
    analyzer = 0
    gc.collect()
    return resultFilename


class timeout:
    """Implements a timeout context manager to work with a with statement."""
    def __init__(self, seconds=1, error_message='Timeout'):
        """

        :param seconds:
        :type seconds: :py:class:`int`
        :param error_message:
        :type error_message: :py:class:`str`
        """
        self.seconds = seconds
        self.error_message = error_message
    def handle_timeout(self, signum, frame):
        raise TimeoutError(self.error_message)
    def __enter__(self):
        signal.signal(signal.SIGALRM, self.handle_timeout)
        signal.alarm(self.seconds)
    def __exit__(self, type, value, traceback):
        signal.alarm(0)
