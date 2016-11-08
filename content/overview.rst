.. _overview:

GZES package overview
=====================

Description
-----------

GZES is a program library for carrying out the equivalent-source processing gravity data. The program library carries out the following functions:

-  The inversion is solved as an optimization problem with the simultaneous goals of (i) minimizing a model objective function and (ii) generating synthetic data that match observations to within a degree of misfit consistent with the statistics of those data.

-  By minimizing the model objective function, distributions of subsurface density contrast are found that are smooth in three dimensions creating a better solution when upward continuing data.

-  The regularization parameter (controlling relative importance of objective function and misfit terms) is determined in either of twoe ways, depending upon how much is known about errors in the measured data.

-  Implementation of parallel computing architecture (OpenMP) allows the user to take full advantage of multi-core processors on a CPU.

- The large size of 3D inversion problems is mitigated by the use of wavelet compression. Parameters controlling the implementation of this compression are available for advanced users. A full sensitivity matrix is also produced to remove the propagation of error into the data space from the lossiness of the wavelet transform.

**NOTE**: The forward modelling of the vertical component of the gravity response to a 3D volume of density contrast can be done after the inversion of a fictitious layer with `gzfor3d`_. 

.. _gzfor3d: http://grav3d.readthedocs.io/en/latest/content/programs/gzfor3d.html

The initial research underlying this program library was funded principally by the mineral industry consortium "Joint and Cooperative Inversion of Geophysical and Geological Data" (1991 - 1997) which was sponsored by NSERC and the following 11 companies: BHP Minerals, CRA Exploration, Cominco Exploration, Falconbridge, Hudson Bay Exploration and Development, INCO Exploration & Technical Services, Kennecott Exploration Company, Newmont Gold Company, Noranda Exploration, Placer Dome, and WMC.

The code itself has been produced and funded by the consortium "Potential fields and software for advanced inversion" (2012-2016) sponsored by Newmont, Teck, Glencore, BHP Billiton, Vale, Computational Geoscience Inc, Cameco, Barrick, Rio Tinto, and Anglo American.

Program library content
-----------------------

Executable programs
^^^^^^^^^^^^^^^^^^^

This package consists of two major programs:

   - GZSENES: Calculates sensitivity matrix
   - GZINVES: Performs inversion of gravity data for an equivalent source layer

Graphical user interfaces
^^^^^^^^^^^^^^^^^^^^^^^^^
GUI-based utilities for these codes include respective viewers for the data and models. They are only available on Windows platforms and can be freely downloaded through the UBC-GIF website:

   - `GM_DATA_VIEWER <http://www.eos.ubc.ca/~rshekhtm/utilities/gm-data-viewer.zip>`__: a utility for viewing raw surface or airborne data (not borehole data), error distributions, and for comparing observed to predicted data directly or as difference maps.
   - `MeshTools3D <http://www.eos.ubc.ca/~rshekhtm/utilities/MeshTools3d.zip>`__: a utility for displaying resulting 3D models as volume renderings. Density contrast volumes can be sliced in any direction, or isosurface renderings can be generated.

Licensing
---------

This equivalent source package is currently only available to the sponsors of "Potential fields and software for advanced inversion" consortium.

.. A **constrained educational version** of the program is available with the `IAG <http://www.flintbox.com/public/project/1605/>`__ package (please visit `UBC-GIF website <http://gif.eos.ubc.ca>`__ for details). The educational version is fully functional so that users can learn how to carry out effective and efficient 3D inversions of magnetic data. **However, RESEARCH OR COMMERCIAL USE IS NOT POSSIBLE because the educational version only allows a limited number of data and model cells**.

.. Licensing for an unconstrained academic version is available - see the `Licensing policy document <http://gif.eos.ubc.ca/software/licenses>`__.

.. **NOTE:** All academic licenses will be **time-limited to one year**. You can re-apply after that time. This ensures that everyone is using the most recent versions of codes.

.. Licensing for commercial use is managed by third party distributors. Details are in the `Licensing policy document <http://gif.eos.ubc.ca/software/licenses>`__.

Installing
----------

There is no automatic installer currently available for the software package. Please follow the following steps in order to use the software:

#. Extract all files provided from the given zip-based archive and place them all together in a new folder such as ``C:\gzes\``

#. Add this directory as new path to your environment variables.

**NOTE**: Do not store anything in the "bin" directory other than executable applications and Graphical User Interface applications (GUIs).


