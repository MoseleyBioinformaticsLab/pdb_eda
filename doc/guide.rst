User Guide
==========

Description
-----------
The :mod:`pdb_eda` package provides a simple Python tool for parsing and analyzing electron density maps data
available from the world wide Protein Data Bank (PDB_).

The :mod:`pdb_eda` package currently provides facilities that can:
    * Parse .ccp4 format file into their object representation.
    * Parse .pdb format file to get information that complimentary to the Bio.PDB module in BioPython_ package.
    * Analyze the electron density maps on atom/residue/chain level and
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

   cd ~/
   git clone https://github.com/MoseleyBioinformaticsLab/pdb_eda.git

Dependencies
~~~~~~~~~~~~

:mod:`pdb_eda` requires the following Python libraries:

   * Biopython_ for creating and analyzing the :mod:`pdb_eda` atom objects.
   * pandas_ for calculating aggregated results.
   * numpy_ and scipy_ for mathmatical calculations.

To install dependencies manually:

.. code:: bash

   pip3 install biopython
   pip3 install pandas
   pip3 install numpy
   pip3 install scipy



Basic usage
-----------
The :mod:`pdb_eda` package can be used in several ways:

   * As a library for accessing and manipulating data in PDB or CCP4 format files.

      * Create the :class:`~pdb_eda.densityAnalysis.fromPDBid` generator function that will generate
        (yield) single :class:`~pdb_eda.densityAnalysis` instance at a time.

      * Process each :class:`~pdb_eda.densityAnalysis` instance:

         * Generate symmetry atoms.
         * Generate red (negative density) or green (positive density) blob lists.
         * Process PDB structures to aggregate cloud.
         * Calculate atom blob list and statistics.

   * As a command-line tool:

      * Calculate statistics of a single PDB structure.
      * Calculate aggregated statistics of multiple PDB structures.





.. _PDB: https://www.wwpdb.org/
.. _BioPython: https://biopython.org/
.. _git: https://git-scm.com/book/en/v2/Getting-Started-Installing-Git/
.. _pandas: http://pandas.pydata.org/
.. _numpy: http://www.numpy.org/
.. _scipy: https://scipy.org/scipylib/index.html