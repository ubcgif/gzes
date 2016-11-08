.. _modelfile:

Model file
==========

This file contains the cell property values (g/cc) of the model and is the most common of the model files. Inversion output models are in this format. The following is the file structure of the model file

.. figure:: ../../images/modelfile.png
    :align: center
    :figwidth: 50%


Each \\(m_{i,j,k}\\) is the property in the \\([i,j,k]^{th}\\) model cell. \\([i, j, k]=[1, 1, 1]\\) is defined as the cell at the top, south-west corner of the model. The total number of lines in this file should equal \\(NN \\times NE \\times NZ\\), where \\(NN\\) is the number of cells in the north direction, \\(NE\\) is the number of cells in the east direction, and \\(NZ\\) is the number of cells in the vertical direction. The model ordering is performed first in the z-direction (top-to-bottom), then in the easting, and finally in the northing.

**NOTE**: Only the cells within the equivalent layer are kept within an inversion. All other cells will be populated with a value of -100. 


