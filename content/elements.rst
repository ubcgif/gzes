Elements of the program GZES
============================

Introduction
------------

The program library consists of the programs:

#. **GZSENES**: calculates sensitivity for the inversion.

#. **GZINVES**: performs gravity equivalent-source processing through inversion.


Each of the above programs requires input files and the specification of parameters in order to run. Although these programs are greatly simplified compared to 3D inversion, we first first present information about general file formats needed.

General files for GZES programs
-------------------------------

There are four general files which are used in GZES. All are in ASCII text format. Input files can have any user-defined name. *Program output files have restricted file names that will be over-written if already in the directory*. Also the filename extensions are not important. Many prefer to use the filename convention so that files are more easily read and edited in the Windows environment. File and file locations may have spaces in the name or path, but it is discouraged. The file name (absolute or relative path) must be 500 characters or less in length. The files contain components of the inversion:

.. toctree::
    :maxdepth: 1

    Modified 3D mesh <files/meshfile>
    Topography <files/topo>
    Observation <files/obs>
    Model (output only) <files/model>

