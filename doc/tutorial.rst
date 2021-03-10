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
    analyzer = densityAnalysis.fromPDBid(pdbid)

The analyzer will only be generated if its .pdb and .ccp4 files exist (valid PDB id),
either locally or can be download on the fly. Or otherwise it will return zero.

Accessing the PDB data
~~~~~~~~~~~~~~~~~~~~~~

The PDB data can be accessed through the **biopdbObj** and **pdbObj**::

    analyzer.biopdbObj
    analyzer.pdbObj

The **biopdbObj** is a Biopython_ data member instance,
and the **pdbObj** is a :class:`pdb_eda.pdbParser.PDBentry` instance that includes some information
that is not available in the Biopython_ instance, such as space group, or rotational matrices.

The information about how to use and access data from the **biopdbObj** instance can be found at Biopython_.

The header information in **pdbObj** can be accessed through *header* attribute as a data member::

   rValue = analyzer.pdbObj.header.rValue
   spaceGroup = analyzer.pdbObj.header.spaceGroup

The available keys include date, method, pdbid, rFree, rValue, resolution, rotationMats, and spaceGroup.
Atom information is optional if  running in *lite* mode.

Accessing the CCP4 data
~~~~~~~~~~~~~~~~~~~~~~~

The CCP4 data can be accessed through the **densityObj** and **diffDensityObj** data members::

    analyzer.densityObj
    analyzer.diffDensityObj

They both contain the header information and the density map from the CCP4 standard map file.
Their header information should be the same, while **densityObj** contains the 2Fo - Fc density map
and **diffDensityObj** contains Fo - Fc density map.
The header information can be accessed through *header* attribute as a data member::

    alpha = analyzer.densityObj.header.alpha
    xlength = analyzer.densityObj.header.xlength

The density map is available in both 1-d and 3-d array::

    oneDmap = analyzer.densityObj.densityArray
    threeDmap = analyzer.densityObj.density

You also have access to several methods that help manipulate the ccp4 data,
for example, to get the point density from a set of given xyz coordinates::

    analyzer.densityObj.getPointDensityFromXyz([10.1, 15.2, 24.4])
    # 1.3517704010009766

A full list of methods can be found at the **API Reference**.

Analyzing the electron density data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There are several methods you can use to perform on the electron density data.
To aggregate the electron density map (2Fo - Fc) by atom, residue, and domain::

    analyzer.aggregateCloud()
    medians = analyzer.medians
    densityElectronRatio = analyzer.densityElectronRatio


To aggregate the difference electron density map (Fo - Fc) into positive (green) and negative (red) blobs::

    greenBlobList = analyzer.greenBlobList
    redBlobList = analyzer.redBlobList

To aggregate the electron density map (2Fo - Fc) into positive (blue) blobs::

    blueBlobList = analyzer.blueBlobList


To acquire a list all nearby symmetry, symmetry-only, or asymmetry atoms::

    symmetryAtoms = analyzer.symmetryAtoms
    symmetryOnlyAtoms = analyzer.symmetryOnlyAtoms
    asymmetryAtoms = analyzer.asymmetryAtoms

To acquire a list all nearby symmetry, symmetry-only, or asymmetry coordinate lists::

    symmetryAtomCoords = analyzer.symmetryAtomCoords
    symmetryOnlyAtomCoords = analyzer.symmetryOnlyAtomCoords
    asymmetryAtomCoords = analyzer.asymmetryAtomCoords

The result is a list of :class:`pdb_eda.densityAnalysis.symAtom` instances.

To calculate the summary statistics of the above positive and negative density blobs with respect to their closest symmetry atom::

    diffMapAtomBlobStatistics = analyzer.calcAtomSpecificBlobStatistics()

For more detailed information, check the **API Reference**.

Using pdb_eda in the command-line interface
-------------------------------------------

Some of the above functions can be accessed from the command line interface::

    Either the "pdb_eda" command or "python3 -m pdb_eda" can be used to run the command line interface.

    > pdb_eda -h

    pdb_eda command-line interface

    Usage:
        pdb_eda -h | --help     for this screen.
        pdb_eda --full-help     help documentation on all modes.
        pdb_eda --version       for the version of pdb_eda.
        pdb_eda single ...      for single structure analysis mode. (Most useful command line mode).
        pdb_eda multiple ...    for multiple structure analysis mode. (Second most useful command line mode).
        pdb_eda contacts ...    for crystal contacts analysis mode.  (Third most useful command line mode).
        pdb_eda generate ...    for generating starting parameters file that then needs to be optimized. (Rarely used mode).
        pdb_eda optimize ...    for parameter optimization mode. (Rarely used mode).

    For help on a specific mode, use the mode option -h or --help.
    For example:
        pdb_eda single --help   for help documentation about single structure analysis mode.



Using single mode to sum significant (> 3 std.dev) deviations in a 3.5 angstrom spherical region around atoms::

   pdb_eda single 3UBK 3ubk.txt difference --atom --radius=3.5 --num-sd=3 --out-format=csv --include-pdbid

Using single mode to sum significant (> 3 std.dev) deviations in a 5 angstrom spherical region around residues::

   pdb_eda single 3UBK 3ubk.txt difference --residue --radius=5 --num-sd=3 --out-format=csv --include-pdbid

Using single mode to return all green difference blobs and their closest symmetry atom::

   pdb_eda single 3UBK 3ubk.green_blobs.txt blob --green --out-format=csv --include-pdbid

Using multiple mode to return summative analysis results for a list of PDB IDs::

   pdb_eda multiple pdbids.txt results/result.txt

Using multiple mode to run single mode with multiprocessing::

   pdb_eda multiple pdbids.txt results/ --single-mode="--atom --radius=3.5 --num-sd=3 --out-format=csv --include-pdbid"

Using multiple mode to check and redownload entry and ccp4 files for a given set of PDB IDs::

   pdb_eda multiple pdbids.txt --reload



.. _PDB: https://www.wwpdb.org/
.. _BioPython: https://biopython.org/
