The pdb_eda Tutorial
====================

The :class:`pdb_eda` package provides classes and other methods for analyzing electron density maps data
available from the worldwide Protein Data Bank (PDB_). It also provides simple command-line interface.


Using pdb_eda as a library
--------------------------

Importing pdb_eda package
~~~~~~~~~~~~~~~~~~~~~~~~~

If the :class:`pdb_eda` is installed, it can be imported::

    import pdb_eda

Constructing densityAnalysis instance
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The densityAnalysis module provides the :func:`~pdb_eda.densityAnalysis.fromPDBid` function that
returns :class:`~pdb_eda.densityAnalysis` instance.
Constructing a :class:`~pdb_eda.densityAnalysis` instance only requires a PDB id::

    pdbid = '1cbs'
    analyser = densityAnalysis.fromPDBid(pdbid)

The analyser will only be generated if its .pdb and .ccp4 files exist (valid PDB id),
either locally or can be download on the fly. Or otherwise it will return zero.

Accessing the PDB data
~~~~~~~~~~~~~~~~~~~~~~

The PDB data can be accessed through the **biopdbObj** and **pdbObj**::

    analyser.biopdbObj
    analyser.pdbObj

The **biopdbObj** is a Biopython_ data member instance,
and the **pdbObj** is a :class:`pdb_eda.pdbParser.PDBentry` instance that includes some information
that is not available in the Biopython_ instance, such as space group, or rotational matrices.

The information about how to use and access data from the **biopdbObj** instance can be found at Biopython_.

The header information in **pdbObj** can be accessed through *header* attribute as a data member::

   rValue = analyser.pdbObj.header.rValue
   spaceGroup = analyser.pdbObj.header.spaceGroup

The available keys include date, method, pdbid, rFree, rValue, resolution, rotationMats, and spaceGroup.
Atom information is optional if  running in *lite* mode.

Accessing the CCP4 data
~~~~~~~~~~~~~~~~~~~~~~~

The CCP4 data can be accessed through the **densityObj** and **diffDensityObj** data members::

    analyser.densityObj
    analyser.diffDensityObj

They both contain the header information and the density map from the CCP4 standard map file.
Their header information should be the same, while **densityObj** contains the 2Fo - Fc density map
and **diffDensityObj** contains Fo - Fc density map.
The header information can be accessed through *header* attribute as a data member::

    alpha = analyser.densityObj.header.alpha
    xlength = analyser.densityObj.header.xlength

The density map is available in both 1-d and 3-d array::

    oneDmap = analyser.densityObj.densityArray
    threeDmap = analyser.densityObj.density

You also have access to several methods that help manipulate the ccp4 data,
for example, to get the point density from a set of given xyz coordinates::

    analyser.densityObj.getPointDensityFromXyz([10.1, 15.2, 24.4])
    # 1.3517704010009766

A full list of methods can be found at the **API Reference**.

Analysing the electron density data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There are several methods you can use to perform on the electron density data.
To aggregate the electron density map (2Fo - Fc) by atom, residue, and chain::

    analyser.aggregateCloud()
    medians = analyser.medians
    chainMedian = analyser.chainMedian


To aggregate the difference electron density map (Fo - Fc) into positive (green) and negative (red) blobs::

    analyser.getBlobList()
    greenBlobList = analyser.greenBlobList
    redBlobList = analyser.redBlobList

To acquire a list all nearby symmetry atoms::

    analyser.calcSymmetryAtoms()
    symmetryAtoms = analyser.symmetryAtoms

The result is a list of :class:`pdb_eda.densityAnalysis.symAtom` instances.

To calculate the summary statistics of the above positive and negative density blobs::

    diffMapStats = analyser.calcAtomBlobDists()

For more detailed information, check the **API Reference**.

Using pdb_eda in the command-line interface
-------------------------------------------

Some of the above functionalities can be accessed as command line interface::

    > python3 -m pdb_eda -h

    pdb_eda command-line interface

    Usage:
        pdb_eda -h | --help
        pdb_eda --version
        pdb_eda single <pdbid> <out-file> [--density-map | --diff-density-map]
        pdb_eda single <pdbid> <out-file> [--radii-param=<paramfile>] [--atom] [--residue] [--chain] [--out-format=<format>]
        pdb_eda single <pdbid> <out-file> [--radii-param=<paramfile>] [--green | --red | --all] [--stats] [--out-format=<format>]
        pdb_eda single <pdbid> <out-file> [--radii-param=<paramfile>] [--symmetry-atoms]
        pdb_eda multiple <pdbid-file> <out-file> [--radii-param=<paramfile>]

    Options:
        -h, --help                      Show this screen.
        --version                       Show version.
        single                          Running single-structure mode
        <pdbid>                         The PDB id
        <out-file>                      Output file name
        multiple                        Running multiple-structure mode
        <pdbid-file>                    File name that contains the pdb ids
        --radii-param=<paramfile>       Radii parameters. [default: conf/optimized_radii_slope_param.json]
        --atom                          Aggregate and print results by atom
        --residue                       Aggregate and print results by residue
        --chain                         Aggregate and print results by chain
        --green                         Calculate and print results of all green blobs (positive difference electron density)
        --red                           Calculate and print results of all red blobs (negative difference electron density)
        --all                           Calculate and print results of both green and red blobs (positive and negative difference electron density)
        --stats                         If set true, return green or red blobs' statistics instead of blob object lists.
        --out-format=<format>           Onput file format, available formats: csv, json [default: json].
        --symmetry-atoms                Calculate and print results of all symmetry atoms. (Only available in jason format)

The single mode will process a single PDB structure and return the desired result file default in json format.
The multiple mode will process multiple PDB structures and return the summary statistics of difference density blobs.

A couple of examples of using the command line interface::

   python3 -m pdb_eda single 3UBK 3ubk.txt --atom --out-format=csv
   python3 -m pdb_eda single 3UBK 3ubk.org.txt --atom --out-format=csv --radii-param='conf/original_radii_slope_param.json'
   python3 -m pdb_eda multiple pdbids.txt results/result.txt


.. _PDB: https://www.wwpdb.org/
.. _BioPython: https://biopython.org/
