# !/usr/bin/python3
"""
pdb.py
    Reads and parses the PDB format files and returns PDB objects.
    Format details of PDB can be found in ftp://ftp.wwpdb.org/pub/pdb/doc/format_descriptions/Format_v33_Letter.pdf

"""
import re
import numpy as np


def readPDBfile(filename, verbose=False):
    """RETURNS DensityMatrix object given the PARAMETER fileName."""
    with open(filename, "r") as fileHandle:
        return parse(fileHandle, verbose)


def parse(handle, mode='lite', verbose=False):
    """RETURNS DensityMatrix object given the PARAMETER file handle."""
    atoms = []
    rotationMats = []
    modelCount = 0
    pdbid = date = method = resolution = rValue = rFree = program = spaceGroup = 0
    for record in handle.readlines():
        if mode == 'lite' and record.startswith('ATOM'):
            break
        elif record.startswith('HEADER'):
            date = record[57: 57+2].strip()
            pdbid = record[62: 62+4].strip()
        elif record.startswith('EXPDTA'):
            method = record[6: 6+30]
            method = method.strip().replace(' ', '_')
        elif record.startswith('REMARK   2 RESOLUTION'):
            match = re.search('RESOLUTION.(.+)ANGSTROMS', record)
            if match:
                resolution = match.group(1).strip()
        elif record.startswith('REMARK   3   R VALUE'):
            match = re.search('^REMARK   3   R VALUE            \(WORKING SET\) : (.+)$', record)
            if match:
                rValue = match.group(1).strip()
        elif record.startswith('REMARK   3   FREE R VALUE'):
            match = re.search('^REMARK   3   FREE R VALUE                     : (.+)$', record)
            if match:
                rFree = match.group(1).strip()
        elif record.startswith('REMARK   3   PROGRAM'):
            match = re.search('^REMARK   3   PROGRAM     : (.+)$', record)
            if match:
                program = match.group(1).strip().replace(' ', '_')
        elif record.startswith('MODEL'):
            modelCount += 1
            if modelCount > 1: break
        elif record.startswith('REMARK 290 SYMMETRY OPERATORS FOR SPACE GROUP:'):
            match = re.search('^REMARK 290 SYMMETRY OPERATORS FOR SPACE GROUP: (.+)$', record)
            if match:
                spaceGroup = match.group(1).strip().replace(' ', '_')
        elif record.startswith('REMARK 290   SMTRY'):
            match = re.search('^REMARK 290   SMTRY(.+)$', record)
            if match:
                items = match.group(1).split()
                if len(rotationMats) < int(items[1]):
                    rotationMats.append(np.zeros((3, 4)))
                rotationMats[int(items[1])-1][int(items[0])-1] = [float(i) for i in items[2:6]]
        elif record.startswith('ATOM') or record.startswith('HETATM'):
            keyValues = {'record': record,
                         'recordType': record[0: 0+6],
                         'serial': record[6: 6+5],
                         'atomName': record[12: 12+4],
                         'alternateLocation': record[16: 16+1],
                         'residueName': record[17: 17+3],
                         'chainID': record[21: 21+1],
                         'residueNumber': record[22: 22+4],
                         'x': record[30: 30+8],
                         'y': record[38: 38+8],
                         'z': record[46: 46+8],
                         'occupancy': record[54: 54+6],
                         'bFactor': record[60: 60+6],
                         'element': record[76: 76+2]}

            keyValues = {key: value.strip() for (key, value) in keyValues.items()}
            atoms.append(Atom(keyValues))


    header = PDBheader(pdbid, date, method, resolution, rValue, rFree, program, spaceGroup, rotationMats)
    return PDBentry(header, atoms)


class PDBentry:
    def __init__(self, header, atoms):
        self.header = header
        self.atoms = atoms


class PDBheader:
    def __init__(self, PDBid, date, method, resolution, rValue, rFree, program, spaceGroup, rotationMats):
        self.pdbid = PDBid
        self.date = date
        self.method = method
        self.resolution = resolution
        self.rValue = rValue
        self.rFree = rFree
        self.program = program
        self.spaceGroup = spaceGroup
        self.rotationMats = rotationMats


class Atom:
    def __init__(self, keyValues):
        self.record = keyValues['record']
        self.recordType = keyValues['recordType']
        self.serial = keyValues['serial']
        self.atomName = keyValues['atomName']
        self.alternateLocation = keyValues['alternateLocation']
        self.residueName = keyValues['residueName']
        self.chainID = keyValues['chainID']
        self.residueNumber = keyValues['residueNumber']
        self.x = keyValues['x']
        self.y = keyValues['y']
        self.z = keyValues['z']
        self.occupancy = keyValues['occupancy']
        self.bFactor = keyValues['bFactor']
        self.element = keyValues['element']




