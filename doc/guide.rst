User Guide
==========

Description
-----------
The :mod:`pdb_eda` package provides a simple Python tool for parsing and analyzing electron density maps data
available from the world wide Protein Data Bank (PDB_).

The :mod:`pdb_eda` package currently provides facilities that can:
    * Parse .ccp4 format file into their object representation.
    * Parse .pdb format file to get information that is complimentary to the Bio.PDB module in BioPython_ package.
    * Analyze the electron density maps at the atom/residue/chain level and
      interpret the electron densities in terms of number of electrons.


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

   * Biopython_ for creating and analyzing the :mod:`pdb_eda` atom objects.
   * numpy_ and scipy_ for mathmatical calculations.

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

    * As a command-line tool:

        * For single-structure mode:
            * Convert electron density map CCP4 files into its equivalent JSON file format.
            * Aggregate electron density map by atom, residue, and chain, and return the results in
              either JSON or csv format.
            * Aggregate difference electron density map into green (positive) or red (negative) blobs,
              and return the object or statistics results in either JSON or csv format.
            * Aggregate difference electron density map for atom and residue specific regions and return
              results in either JSON or csv format.
            * Return traditional quality metrics and statistics for atoms and residues.

        * For multiple-structure mode:
            * Analyze and return cumulative statistics for a given list of PDB IDs.
            * Filter list of PDB IDs by cumulative statistic criteria.
            * Check and redownload problematic PDB entries.
            * Run single structure mode with multicore processing.
            * Run crystal contacts mode with multicore processing.

        * For crystal contacts mode:
            * Analyze and return atoms with crystal contacts.
            * This mode requires pymol to be installed.

        * For parameter generation mode: (rarely used mode)
            * Downloads PDB chemical component list and extracts information to create atom type parameters.
            * Analyzes list of PDB IDs for specific atom types.
            * Generates atom type parameter file and list of PDB IDs for their optimization.

        * For parameter optimization mode: (rarely used mode)
            * Optimizes atom type radii and b-factor density correction slopes using a given list of PDB IDs.

.. _PDB: https://www.wwpdb.org/
.. _BioPython: https://biopython.org/
.. _git: https://git-scm.com/book/en/v2/Getting-Started-Installing-Git/
.. _numpy: http://www.numpy.org/
.. _scipy: https://scipy.org/scipylib/index.html