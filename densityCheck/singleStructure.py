import sys
import math
import numpy as np
from scipy import stats
import multiprocessing
import datetime

from . import densityAnalysis

## Final set of radii derived from optimization
radiiDefault = {'C_single': 0.84, 'C_double': 0.67, 'C_intermediate': 0.72, 'C_single_bb': 0.72, 'C_double_bb': 0.61, 
                'O_single': 0.80, 'O_double': 0.77, 'O_intermediate': 0.82, 'O_double_bb': 0.71, 
                'N_single': 0.95, 'N_intermediate': 0.77, 'N_single_bb': 0.7, 
                'S_single': 0.75}

slopesDefault = {'C_single': -0.33373359301010996, 'C_double': -0.79742538209281033, 'C_intermediate': -0.46605044936397311, 'C_single_bb': -0.44626029983492604, 'C_double_bb': -0.53313491315106321, 
                'O_single': -0.61612691527685515, 'O_double': -0.69081386892018048, 'O_intermediate': -0.64900152439990022, 'O_double_bb': -0.59780395867126568, 
                'N_single': -0.50069328952200343, 'N_intermediate': -0.5791543104296577, 'N_single_bb': -0.4905830761073553, 
                'S_single': -0.82192528851026947}


def main(pdbidfile, pdbid, resultname):
    '''
    Mode: 0 - use optimized radii
          1 - use original radii
    '''
    analyser = densityAnalysis.fromPDBid(pdbid)
    analyser.aggregateCloud(radiiDefault, slopesDefault, atomList=True)

    fileHandle = open(resultname, 'w')
    for item in analyser.atomList:
        print(', '.join(map(str, item + [analyser.chainMedian])), file=fileHandle)

if __name__ == '__main__':
    if len(sys.argv) == 3:
        _, pdbid, resultname = sys.argv
        main(filename, pdbid, resultname)
    else:
        print("Wrong number of arguments!")

