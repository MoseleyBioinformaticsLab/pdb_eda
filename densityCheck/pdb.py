# !/usr/bin/python3
"""
pdb.py
    Reads and parses the PDB format files and returns PDB objects.
    Format details of PDB can be found in ftp://ftp.wwpdb.org/pub/pdb/doc/format_descriptions/Format_v33_Letter.pdf

"""
import gzip
import re


def read(filename, verbose=False):
    """RETURNS DensityMatrix object given the PARAMETER fileName."""
    with gzip.open(filename, "rb") as fileHandle:
        return parse(fileHandle, verbose)


def parse(handle, verbose=False):
    """RETURNS DensityMatrix object given the PARAMETER file handle."""
    atoms = []
    modelCount = 0
    for record in handle:
        if record.startswith('ATOM') or record.startswith('HETATM'):
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
                         'element': record[76: 76+2],
                         'PDBid': PDBid,
                         'method': method,
                         'date': date,
                         'resolution': resolution,
                         'rValue': rValue,
                         'rFree': rFree}

            keyValues = {key: value.strip() for (key, value) in keyValues.items()}
            atoms.append(Atom(keyValues))
        elif record.startswith('HEADER'):
            PDBid = record[62: 62+4]
            date = record[57, 57+2]
        elif record.startswith('EXPDATA'):
            method = record[6: 6+30]
            method = method.strip()
            method = method.replace(' ', '_')
        elif record.startswith('REMARK   2 RESOLUTION'):
            match = re.search('RESOLUTION.(.+)ANGSTROMS', record)
            if match:
                resolution = match.group(1)
        elif record.startswith('REMARK   3   R VALUE'):
            match = re.search('^REMARK   3   R VALUE            \(WORKING SET\) : (.+)$', record)
            if match:
                rValue = match.group(1)
        elif record.startswith('REMARK   3   FREE R VALUE'):
            match = re.search('^REMARK   3   FREE R VALUE                     : (.+)$', record)
            if match:
                rFree = match.group(1)
        elif record.startswith('MODEL'):
            modelCount += 1
            if modelCount >1: break




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
        self.file = keyValues['file']
        self.PDBid = keyValues['PDBid']
        self.method = keyValues['method']
        self.date = keyValues['date']
        self.resolution = keyValues['resolution']
        self.rValue = keyValues['rValue']
        self.rFree = keyValues['rFree']
        self.solvent = keyValues['solvent']


