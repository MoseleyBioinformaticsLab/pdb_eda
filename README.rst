pdb-eda
==========

.. image:: https://raw.githubusercontent.com/MoseleyBioinformaticsLab/pdb_eda/master/doc/_static/images/pdb_eda_logo.png
   :width: 50%
   :align: center
   :target: https://pdb-eda.readthedocs.io/

Description
-----------
The `pdb_eda` package provides a simple Python tool for parsing and analyzing electron density maps data
available from the world wide Protein Data Bank (PDB_).

The `pdb_eda` package currently provides facilities that can:
    * Parse .ccp4 format file into their object representation.
    * Parse .pdb format file to get information that complimentary to the Bio.PDB module in BioPython_ package.
    * Analyze the electron density maps on atom/residue/chain level and
      interpret the electron densities in terms of number of electrons.

Full API documentation, user guide, and tutorial can be found on readthedocs_

Installation
------------
`pdb_eda` runs under Python 3.4+ and is available through python3-pip.
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

`pdb_eda` requires the following Python libraries:

   * Biopython_ for creating and analyzing the `pdb_eda` atom objects.
   * pandas_ for calculating aggregated results.
   * numpy_ and scipy_ for mathmatical calculations.
   * docopt_ for better command line interface.
   * jsonpickle_ for formatted and reusable output.

To install dependencies manually:

.. code:: bash

   pip3 install biopython
   pip3 install pandas
   pip3 install numpy
   pip3 install scipy
   pip3 install docopt
   pip3 install jsonpickle


Basic usage
-----------
The `pdb_eda` package can be used in several ways:

   * As a library for accessing and manipulating data in PDB or CCP4 format files.

      * Create the `pdb_eda.densityAnalysis.fromPDBid` generator function that will generate
        (yield) single `pdb_eda.densityAnalysis` instance at a time.

      * Process each `pdb_eda.densityAnalysis` instance:

         * Generate symmetry atoms.
         * Generate red (negative density) or green (positive density) blob lists.
         * Process PDB structures to aggregate cloud.
         * Calculate atom blob list and statistics.

   * As a command-line tool:

      * Calculate statistics of a single PDB structure.
      * Calculate aggregated statistics of multiple PDB structures.




.. _readthedocs: https://pdb-eda.readthedocs.io/en/latest/
.. _PDB: https://www.wwpdb.org/
.. _BioPython: https://biopython.org/
.. _git: https://git-scm.com/book/en/v2/Getting-Started-Installing-Git/
.. _pandas: http://pandas.pydata.org/
.. _numpy: http://www.numpy.org/
.. _scipy: https://scipy.org/scipylib/index.html
.. _docopt: http://docopt.org/
.. _jsonpickle: https://github.com/jsonpickle/jsonpickle
