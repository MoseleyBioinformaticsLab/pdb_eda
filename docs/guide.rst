User Guide
==========

Description
-----------
The :mod:`pdb_eda` package provides a simple Python tool for parsing and analyzing electron density maps data
available from the world wide Protein Data Bank (PDB_).

The :mod:`pdb_eda` package currently provides facilities that can:
    * Parse .ccp4 format file into their object representation.
    * Parse .pdb format file to get information that is complimentary to the Bio.PDB module in BioPython_ package.
    * Analyze the electron density maps at the atom/residue/domain levels and
      interpret the electron densities in terms of number of electrons by estimating a density-electron ratio.

Citation
--------
Please cite the following papers when using pdb_eda:

Sen Yao and Hunter N.B. Moseley. "A chemical interpretation of protein electron density maps in the worldwide protein data bank" PLOS One 15, e0236894 (2020).
https://doi.org/10.1371/journal.pone.0236894

Sen Yao and Hunter N.B. Moseley. "Finding high-quality metal ion-centric regions across the worldwide Protein Data Bank" Molecules 24, 3179 (2019).
https://doi.org/10.3390/molecules24173179

Installation
------------
:mod:`pdb_eda` runs under Python 3.4+ and is available through python3-pip.
Install via pip or clone the git repo and install the following dependencies and you are ready to go!

Install on Linux, Mac OS X
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: bash

   python3 -m pip install pdb_eda

GitHub Package installation
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Make sure you have git_ installed:

.. code:: bash

   git clone https://github.com/MoseleyBioinformaticsLab/pdb_eda.git

Dependencies
~~~~~~~~~~~~

:mod:`pdb_eda` requires the following Python libraries:

   * Biopython_ for creating and analyzing the `pdb_eda` atom objects.
   * Cython_ for cythonizing low-level utility functions to improve computational performance.
      * Requires gcc to be installed for the cythonization process.
   * numpy_ and scipy_ for mathmatical calculations.
   * docopt_ for better command line interface.
   * jsonpickle_ for formatted and reusable output.
   * PyCifRW_ for reading Cif formatted files.
      * Requires gcc to be installed for compiling components of the package.
   * pymol_ for calculating crystal contacts. (This package is not required, except for this functionality).

To install dependencies manually:

.. code:: bash

   pip3 install biopython
   pip3 install cython
   pip3 install numpy
   pip3 install scipy
   pip3 install docopt
   pip3 install jsonpickle
   pip3 install PyCifRW


Basic usage
-----------
The :mod:`pdb_eda` package can be used in several ways:

    * As a library for accessing and manipulating data in PDB and CCP4 format files.

        * Use the :class:`~pdb_eda.densityAnalysis.fromPDBid` generator function that will generate
          (yield) a single :class:`~pdb_eda.densityAnalysis` instance at a time.
        * Process each :class:`~pdb_eda.densityAnalysis` instance:
        * Generate symmetry atoms.
        * Generate red (negative density) or green (positive density) blob lists.
        * Process PDB structures to aggregate cloud.
        * Calculate atom blob list and statistics.
        * Calculate atom regional discrepancies and statistics.
        * Calculate residue regional discrepancies and statistics.

    * As a command-line tool using the pdb_eda command (or "python3 -m pdb_eda").

        * The command-line interface has multiple modes.

        * single - single-structure mode:
            * Convert electron density map CCP4 files into its equivalent JSON file format.
            * Aggregate electron density map by atom, residue, and domain, and return the results in
              either JSON or csv format.
            * Aggregate difference electron density map into green (positive) or red (negative) blobs,
              and return the object or statistics results in either JSON or csv format.
            * Aggregate difference electron density map for atom and residue specific regions and return
              results in either JSON or csv format.
            * Return traditional quality metrics and statistics for atoms and residues.

        * multiple - multiple-structure mode:
            * Analyze and return cumulative statistics for a given list of PDB IDs.
            * Filter list of PDB IDs by cumulative statistic criteria.
            * Check and redownload problematic PDB entries.
            * Run single structure mode with multicore processing.
            * Run crystal contacts mode with multicore processing.

        * contacts - crystal contacts mode:
            * Analyze and return atoms with crystal contacts.
            * This mode requires pymol to be installed.

        * generate - parameter generation mode: (rarely used mode)
            * Downloads PDB chemical component list and extracts information to create atom type parameters.
            * Analyzes list of PDB IDs for specific atom types.
            * Generates atom type parameter file and list of PDB IDs for their optimization.

        * optimize - parameter optimization mode: (rarely used mode)
            * Optimizes atom type radii and b-factor density correction slopes using a given list of PDB IDs.


CHANGELOG
---------
Since version 1.0.1, over 2200 lines of additional code has been written and most of the code base has been revised and refactored.
Computationally intensive parts of the code have been cythonized to improve execution performance.
Many variables and functions have been renamed to greatly improve readability and understanding of the code base, API, and CLI.

The application programming interface (API) has been greatly expanded and much of the functionality streamlined.

The command line interface has been greatly expanded and now includes single, multiple, contacts, generate, and optimize modes.

Optimize mode has a new penalty function being optimized that both minimizes differences in density-electron ratio estimates and
maximizes electron cloud aggregation.  The optimization is also roughly 10-fold faster than the previous generation of algorithm.

The atom types have been systematically generated from the wwPDB master chemical components file.
Both amino acid and nucleic acid type parameters have been optimized.
So both protein and nucleic acid PDB entries can be analyzed now.


License
-------
A modified Clear BSD License

Copyright (c) 2019, Sen Yao, Hunter N.B. Moseley
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted (subject to the limitations in the disclaimer
below) provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

* Neither the name of the copyright holder nor the names of its contributors may be used
  to endorse or promote products derived from this software without specific
  prior written permission.

* If the source code is used in a published work, then proper citation of the source
  code must be included with the published work.

NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE GRANTED BY THIS
LICENSE. THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.

.. _readthedocs: https://pdb-eda.readthedocs.io/en/latest/
.. _PDB: https://www.wwpdb.org/
.. _BioPython: https://biopython.org/
.. _Cython: https://cython.readthedocs.io/en/latest/index.html
.. _git: https://git-scm.com/book/en/v2/Getting-Started-Installing-Git/
.. _numpy: http://www.numpy.org/
.. _scipy: https://scipy.org/scipylib/index.html
.. _docopt: http://docopt.org/
.. _jsonpickle: https://github.com/jsonpickle/jsonpickle
.. _PyCifRW: https://pypi.org/project/PyCifRW/4.3/
.. _pymol: https://pymol.org/2/
