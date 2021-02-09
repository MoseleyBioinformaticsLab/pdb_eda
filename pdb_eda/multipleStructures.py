"""
pdb_eda multiple structure analysis mode command-line interface
  Analyzes multiple pdb entries.

Usage:
    pdb_eda multiple -h | --help
    pdb_eda multiple <pdbid-file> <out-result-file> [--params=<params-file>] [--out-format=<format>] [--testing] [--time-out=<seconds>]
    pdb_eda multiple <in-result-file> <out-pdbid-file> --filter [--out-format=<format>] [--max-resolution=<max-resolution>] [--min-resolution=<min-resolution>] [--min-atoms=<min-atoms>] [--min-residues=<min-residues>] [--min-elements=<min-elements>]
    pdb_eda multiple <pdbid-file> <out-dir> --single-mode=<quoted-single-mode-options> [--time-out=<seconds>] [--testing]
    pdb_eda multiple <pdbid-file> <out-dir> --contacts-mode=<quoted-contacts-mode-options> [--time-out=<seconds>] [--testing]

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
    --single-mode=<quoted-single-mode-options>      Run single structure analysis mode on a set of PDB entries.
    --contacts-mode=<quoted-contacts-mode-options>  Run (crystal) contacts analysis mode on a set of PDB entries.
    --filter                                        Filter pdbids based on criteria.
    --max-resolution=<max-resolution>               Maximum x-ray crystallographic resolution to allow. [default: 3.5]
    --min-resolution=<min-resolution>               Minimum x-ray crystallographic resolution to allow. [default: 0]
    --min-atoms=<min-atoms>                         Minimum number of atoms analyzed (numAtomsAnalyzed) required. [default: 1000]
    --min-residues=<min-residues>                   Minimum number of residues required of the given residue types separated by commas. [default: 0]
    --min-elements=<min-elements>                   Minimum number of atoms required of the given elements separated by commas. [default: 0]
"""

import docopt
import os
import time
import sys

import json
import csv
import multiprocessing
import tempfile
import signal
import collections

from . import densityAnalysis
from . import __version__
from . import singleStructure
from . import crystalContacts

defaultParamsFilepath = os.path.join(os.path.dirname(__file__), 'conf/optimized_params.json')

globalParams = None
globalArgs = {}

def main():
    global globalArgs
    globalArgs = docopt.docopt(__doc__, version=__version__)
    if globalArgs["--help"]:
        print(__doc__)
        exit(0)
    globalArgs["--time-out"] = int(globalArgs["--time-out"])

    paramsFilepath = globalArgs["--params"] if globalArgs["--params"] else defaultParamsFilepath
    try:
        with open(paramsFilepath, 'r') as paramsFile:
            global globalParams
            globalParams = json.load(paramsFile)
    except:
        sys.exit(str("Error: params file \"") + paramsFilepath + "\" does not exist or is not parsable.")

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
            sys.exit(str("Error: PDB results file \"") + globalArgs["<in-result-file>"] + "\" does not exist or is not parsable.")

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
    else:
        try:
            pdbids = []
            with open(globalArgs["<pdbid-file>"], "r") if globalArgs["<pdbid-file>"] != "-" else sys.stdin as textFile:
                for pdbid in textFile:
                    pdbids.append(pdbid[0:4])
        except:
            sys.exit(str("Error: PDB IDs file \"") + globalArgs["<pdbid-file>"] + "\" does not exist or is not parsable.")

        if globalArgs["--single-mode"]:
            processFunction = singleModeFunction
            if not os.path.isdir(globalArgs["<out-dir>"]):
                if not os.path.isfile(globalArgs["<out-dir>"]):
                    os.mkdir(globalArgs["<out-dir>"])
                else:
                    sys.exit(str("Error: Output directory \"") + globalArgs["<out-dir>"] + "\" is a file.")
        elif globalArgs["--contacts-mode"]:
            processFunction = contactsModeFunction
            if not os.path.isdir(globalArgs["<out-dir>"]):
                if not os.path.isfile(globalArgs["<out-dir>"]):
                    os.mkdir(globalArgs["<out-dir>"])
                else:
                    sys.exit(str("Error: Output directory \"") + globalArgs["<out-dir>"] + "\" is a file.")
        else:
            processFunction = multipleModeFunction

        if globalArgs["--testing"]:
            results = [ processFunction(pdbid) for pdbid in pdbids ]
        else:
            with multiprocessing.Pool() as pool:
                results = pool.map(processFunction, pdbids)

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
                statsHeaders = ['density_electron_ratio', 'voxel_volume', 'f000', 'chain_num_voxel', 'chain_total_electrons', 'density_mean', 'diff_density_mean', 'resolution', 'space_group', 'num_atoms_analyzed', 'num_residues_analyzed', 'num_chains_analyzed']
                with open(globalArgs['<out-result-file>'], "w", newline='') if globalArgs["<out-result-file>"] != "-" else sys.stdout as csvFile:
                    writer = csv.writer(csvFile)
                    writer.writerow(['pdbid'] + statsHeaders + sorted(globalParams["radii"]))
                    for result in fullResults.values():
                        stats = [result["stats"][header] for header in statsHeaders]
                        diffs = [result["diffs"][atomType] for atomType in sorted(globalParams["radii"])]
                        writer.writerow([result['pdbid']] + stats + diffs)
            else:
                with open(globalArgs['<out-result-file>'], "w") if globalArgs["<out-result-file>"] != "-" else sys.stdout as jsonFile:
                    print(json.dumps(fullResults, indent=2, sort_keys=True), file=jsonFile)


def singleModeFunction(pdbid):
    """
    Execute pdb_eda single mode on a given pdbid.

    :param str pdbid:
    :return: 0
    """
    singleModeCommandLine = "pdb_eda single " + pdbid + " " + globalArgs["<out-dir>"] + "/" + pdbid + ".result " + globalArgs["--single-mode"]
    sys.argv = singleModeCommandLine.split()
    if globalArgs["--time-out"]:
        try:
            with timeout(seconds=globalArgs["--time-out"]):
                singleStructure.main()
        except:
            pass
    elif globalArgs["--testing"]:
        singleStructure.main()
    else:
        try:
            singleStructure.main()
        except:
            pass

    return 0


def contactsModeFunction(pdbid):
    """
    Execute pdb_eda contacts mode on a given pdbid.

    :param str pdbid:
    :return: 0
    """
    contactsModeCommandLine = "pdb_eda contacts " + pdbid + " " + globalArgs["<out-dir>"] + "/" + pdbid + ".result " + globalArgs["--contacts-mode"]
    sys.argv = contactsModeCommandLine.split()
    if globalArgs["--time-out"]:
        try:
            with timeout(seconds=globalArgs["--time-out"]):
                crystalContacts.main()
        except:
            pass
    elif globalArgs["--testing"]:
        crystalContacts.main()
    else:
        try:
            crystalContacts.main()
        except:
            pass

    return 0


def multipleModeFunction(pdbid):
    """
    Analyze a single pdbid.

    :param str pdbid:
    :return: resultFilename or 0
    """
    if globalArgs["--time-out"]:
        try:
            with timeout(seconds=globalArgs["--time-out"]):
                return analyzePDBID(pdbid)
        except:
            return 0
    else:
        return analyzePDBID(pdbid)


def analyzePDBID(pdbid):
    """Process function to analyze a single pdb entry.

    :param :py:class:`str` pdbid: pdbid for entry to download and analyze.
    :return: resultFilename
    :rtype: :py:class:`str`
    """
    startTime = time.process_time()

    analyzer = densityAnalysis.fromPDBid(pdbid)

    if not analyzer:
        return 0

    analyzer.aggregateCloud(globalParams, atomL=True, residueL=True, chainL=True)
    analyzer.estimateF000()
    if not analyzer.densityElectronRatio:
        return 0

    diffs = { atomType:((analyzer.medians['corrected_density_electron_ratio'][atomType] - analyzer.densityElectronRatio) / analyzer.densityElectronRatio)
              if atomType in analyzer.medians['corrected_density_electron_ratio'] else 0 for atomType in sorted(globalParams["radii"]) }

    stats = { 'density_electron_ratio' : analyzer.densityElectronRatio, 'voxel_volume' : analyzer.densityObj.header.unitVolume, 'f000' : analyzer.f000, 'chain_num_voxel' : analyzer.chainNvoxel,
        'chain_total_electrons' : analyzer.chainTotalE, 'density_mean' : analyzer.densityObj.header.densityMean, 'diff_density_mean' : analyzer.diffDensityObj.header.densityMean,
        'resolution' : analyzer.pdbObj.header.resolution, 'space_group' : analyzer.pdbObj.header.spaceGroup, 'num_atoms_analyzed' : len(analyzer.atomList),
        'num_residues_analyzed' : len(analyzer.residueList), 'num_chains_analyzed' : len(analyzer.chainList)  }

    properties = { property : value for (property,value) in analyzer.biopdbObj.header.items() }
    properties['residue_counts'] = dict(collections.Counter(residue.resname for residue in analyzer.biopdbObj.get_residues()))
    properties['element_counts'] = dict(collections.Counter(atom.element for atom in analyzer.biopdbObj.get_atoms()))

    elapsedTime = time.process_time() - startTime
    resultFilename = createTempJSONFile({ "pdbid" : analyzer.pdbid, "diffs" : diffs, "stats" : stats, "execution_time" : elapsedTime, "properties" : properties }, "tempResults_")
    return resultFilename


def createTempJSONFile(data, filenamePrefix):
    """Creates a temporary JSON file and returns its filename.

    :param data:  data to save into the JSON file.
    :param :py:class:`str` filenamePrefix: temporary filename prefix.
    :return: filename
    :rtype: :py:class:`str`
    """
    dirname = os.getcwd()
    filename = 0
    with tempfile.NamedTemporaryFile(mode='w', buffering=1, dir=dirname, prefix=filenamePrefix, delete=False) as tempFile:
        json.dump(data,tempFile)
        filename = tempFile.name
    return filename


class timeout:
    """
    Implements a timeout context manager to work with a with statement.
    """
    def __init__(self, seconds=1, error_message='Timeout'):
        self.seconds = seconds
        self.error_message = error_message
    def handle_timeout(self, signum, frame):
        raise TimeoutError(self.error_message)
    def __enter__(self):
        signal.signal(signal.SIGALRM, self.handle_timeout)
        signal.alarm(self.seconds)
    def __exit__(self, type, value, traceback):
        signal.alarm(0)
