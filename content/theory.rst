.. _theory:

Background theory
=================

Introduction
------------

.. _grav3d: http://grav3d.readthedocs.io

The GZES suite of algorithms is based on the grav3d_ suit, developed at the UBC Geophysical Inversion Facility, is used to invert gravimetric responses over a three dimensional distribution of density contrast, or anomalous density. This manual is designed so that geophysicists who are familiar with the gravity experiment, but who are not necessarily versed in the details of inverse theory, can use the codes and invert their data. In the following, we describe the basics of the algorithm, but readers are referred to :cite:`LiOldenburg10` for an in-depth discussion of various aspects of the algorithm. Note that an understanding of these components is necessary for the user to have a global view of the algorithms and to use the program library to its fullest extent. 

A gravity experiment involves measuring the vertical components of the gravity field produced by anomalous (either excess or deficient) mass beneath the surface. A distribution of anomalous mass, characterized by anomalous density :math:`\rho(x, y, z)`, produces its own gravity field, :math:`\mathbf{g}_s`, which is superimposed on the ambient gravity field. By measuring the resultant field and removing the ambient field from the measurements through numerical processing, one obtains the field due to the anomalous mass.

The vertical component of the gravity field produced by the density :math:`\rho(x, y, z)` is given by

.. _gzfield_:
.. math:: 
    g_z(\mathbf{r}_o)= \gamma\int\limits_V \rho(\mathbf{r})\frac{z-z_o}{\left | \mathbf{r}-\mathbf{r}_o \right |^3} dv,
    :label: gzfield

where :math:`\mathbf{r}_o = (x_o,y_o,z_o)` is the vector denoting the observation location and :math:`\mathbf{r} = (x,y,z)` is the source location. The volume of the anomalous mass is :math:`V` and :math:`\gamma` is the gravitational constant. Here we have adopted a Cartesian coordinate system having its origin on the earth's surface and the :math:`z-`\ axis pointing vertically downward. In the following, we outline the basics of the forward and inverse procedures used by the GRAV3D program library.


Forward Modelling
-----------------

Forward modelling of gravity data is a linear problem and can be carried out by performing the integration in equation :eq:`gzfield`. We divide the region of interest into a set of 3D prismatic cells by using a 3D orthogonal mesh with a single layer and assume a constant density contrast within each cell. If topography is present, the layer is draped on topography. Also, We discretize the density contrast model in this manner since it is best suited for our inversion methodology. Given such a discretization, the gravity field at the :math:`i^{th}` location can be written as:

.. _fwdata_:

.. math:: 
     \begin{aligned}
     d_i\equiv g_z(\mathbf{r}_{oi}) \\
     = \overset{M}{\underset{j=1}{\sum}}\rho_j\left \{ \gamma\int\limits_{\Delta V_j}\frac{z-z_o}{\left | \mathbf{r}-\mathbf{r}_{oi} \right |^3}dv \right \}, \\
     \equiv \overset{M}{\underset{j=0}{\sum}}\rho_j G_{ij}.
     \end{aligned}
     :label: fwdata


In equation :eq:`fwdata`, :math:`\rho_j` and :math:`\Delta V_j` are the anomalous density and volume of the :math:`j^{th}` cell, :math:`d_i` is introduced as a generic symbol for the :math:`i^{th}` datum, and :math:`G_{ij}`, defined by the expression in brackets, quantifies the contribution of the :math:`j^{th}` cell to the :math:`i^{th}` datum. The solution for the integral in equation :eq:`fwdata` can be found in :cite:`Nagy66` and we have adopted the solution by :cite:`Haaz53` here.

Inversion methodology
---------------------

Let the set of extracted anomaly data be :math:`\mathbf{d} = (d_1,d_2,...,d_N)^T` and the density contrast of cells in the model be :math:`\mathbf{\rho} = (\rho_1,\rho_2,...,\rho_M)^T`. The two are related by the sensitivity matrix

.. _sens_:
.. math::
     \mathbf{d}=\mathbf{G}\mathbf{\rho}.
     :label: sens

The matrix has elements :math:`g_{ij}` which quantify the contribution to the :math:`i^{th}` datum due to a unit density in the :math:`j^{th}` cell. The program performs the calculation of the sensitivity matrix, which is to be used by the subsequent inversion. The sensitivity matrix provides the forward mapping from the model to the data during the entire inverse process. We will discuss its efficient representation via the wavelet transform in a separate section. 

For the inversion, the first question that arises concerns definition of the "model". We choose density contrast, :math:`\rho`, as the model for since the anomalous field is directly proportional to the density contrast. The inverse problem is formulated as an optimization problem where a global objective function, :math:`\phi`, is minimized subject to the constraints in equation :eq:`sens`. The global objective functions consists of two components: a model objective function, :math:`\phi_m`, and a data misfit function, :math:`\phi_d`, such that

.. _globphi_:
.. math::
    \begin{aligned}
    \min \phi = \phi_d+\beta\phi_m \\
    \mbox{s. t. } \rho^l\leq \rho \leq \rho^u, \nonumber
    \end{aligned}
    :label: globphi

where :math:`\beta` is a trade off parameter that controls the relative importance of the model smoothness through the model objective function and data misfit function. When the standard deviations of data errors are known, the acceptable misfit is given by the expected value :math:`\phi_d` and we will search for the value of :math:`\beta` via an L-curve criterion :cite:`Hansen00` that produces the expected misfit. Otherwise, a user-defined value is used. Bound are imposed through the projected gradient method so that the recovered model lies between imposed lower (:math:`\rho^l`) and upper (:math:`\rho^u`) bounds.

We next discuss the construction of a model objective function which, when minimized, produces a model that is geophysically interpretable. However, here we stress that the functional is to create a smooth model since the layer of sources being solved for is completely fictional and a means to create a model that re-produces the data. In general, the objective function gives the flexibility to incorporate as little or as much information as possible. At the very minimum, this function drives the solution towards a reference model :math:`\rho_o` and requires that the model be relatively smooth in the three spatial directions. Here we adopt a right handed Cartesian coordinate system with positive north and positive down. Let the model objective function be

.. _mof:
.. math::
    \phi_m(\rho) &=& \alpha_s\int\limits_V \left\{\rho(\mathbf{r}) \right\}^2dv + \alpha_x\int\limits_V \left\{\frac{\partial \rho(\mathbf{r})}{\partial x}\right\}^2dv + \alpha_y\int\limits_V \left\{\frac{\partial \rho(\mathbf{r})}{\partial y}\right\}^2dv
    :label: mof

where :math:`\alpha_s`, :math:`\alpha_x`, and :math:`\alpha_y` are coefficients, which affect the relative importance of different components in the objective function. Numerically, the model objective function in equation :eq:`mof` is discretized onto the mesh defining the density contrast model using a finite difference approximation. This yields:

.. _modobjdiscr_:
.. math::
    \begin{aligned}
    \phi_m(\mathbf{\rho}) = \mathbf{\rho}^T(\alpha_s \mathbf{W}_s^T\mathbf{W}_s + \alpha_x \mathbf{W}_x^T\mathbf{W}_x+\alpha_y \mathbf{W}_y^T\mathbf{W}_y)\mathbf{\rho}, \nonumber\\
    \equiv \mathbf{\rho}^T\mathbf{W}_m^T\mathbf{W}_m\mathbf{\rho}, \nonumber\\
    =\left \| \mathbf{W}_m \mathbf{\rho} \right \|^2,\end{aligned}
    :label: modobjdiscr


where :math:`\mathbf{\rho}` is an :math:`M`-length vector representing the recovered model. The next step in setting up the inversion is to define a misfit measure. Here we use the :math:`l_2`-norm measure

.. _phid_:
.. math::
    \phi_d = \left\| \mathbf{W}_d(\mathbf{G}\mathbf{\rho}-\mathbf{d})\right\|^2.
    :label: phid

For the work here, we assume that the contaminating noise on the data is independent and Gaussian with zero mean. Specifying :math:`\mathbf{W}_d` to be a diagonal matrix whose :math:`i^{th}` element is :math:`1/\sigma_i`, where :math:`\sigma_i` is the standard deviation of the :math:`i^{th}` datum makes :math:`\phi_d` a chi-squared distribution with :math:`N` degrees of freedom. The optimal data misfit for data contaminated with independent, Gaussian noise has an expected value of :math:`E[\chi^2]=N`, providing a target misfit for the inversion. We now have the components to solve the inversion as defined in equation :eq:`globphi`.

To solve the optimization problem when constraints are imposed we use the projected gradients method :cite:`CalamaiMore87,Vogel02`. This technique forces the gradient in the Krylov sub-space minimization (in other words a step during the conjugate gradient process) to zero if the proposed step would make a model parameter exceed the bound constraints. The result is a model that reaches the bounds, but does not exceed them. This method is computationally faster than the log-barrier method because (1) model parameters on the bounds are neglected for the next iteration and (2) the log-barrier method requires the calculation of a barrier term. Previous versions of potential-field codes from UBC-GIF used the logarithmic barrier method :cite:`Wright97,NocedalWright99`.

.. _waveletSection:

Wavelet Compression of Sensitivity Matrix
-----------------------------------------

The major obstacle to the solution of a large-scale gravity processing problem is the CPU time required for the application of the sensitivity matrix to model vectors. The GZES program library overcomes these difficulties by forming a sparse representation of the sensitivity matrix using a wavelet transform based on compactly supported, orthonormal wavelets. For more details, the users are referred to :cite:`LiOldenburg03,LiOldenburg10`. In the following, we give a brief description of the method necessary for the use of the GZES library. 

Each row of the sensitivity matrix in a 3D gravity inversion can be treated as a 2D image and a 2D wavelet transform can be applied to it. By the properties of the wavelet transform, most transform coefficients are nearly or identically zero. When coefficients of small magnitudes are discarded (the process of thresholding), the remaining coefficients still contain much of the necessary information to reconstruct the sensitivity accurately. These retained coefficients form a sparse representation of the sensitivity in the wavelet domain. The need to store only these large coefficients means that the memory requirement is reduced. Further, the multiplication of the sensitivity with a vector can be carried out by a sparse multiplication in the wavelet domain. This greatly reduces the CPU time. Since the matrix-vector multiplication constitutes the core computation of the inversion, the CPU time for the inverse solution is reduced accordingly. The use of this approach increases the size of solvable problems by nearly two orders of magnitude. However, the model subtleties from the wavelet compression will have an impact on the forward modelled data (typically these are not worried about when 3D inversion is performed as the minute differences in the models won't affect intpertation). Therefore, the full sensitivity matrix is also written to file and can be used at the end of the inversion allowing the inversion code to solve the problem almost the entire way in the wavelet domain and then finish with the accuracy needed to create an accurate representation of the data.

Let :math:`\mathbf{G}` be the sensitivity matrix and :math:`\mathcal{W}` be the symbolic matrix-representation of the 3D wavelet transform. Then applying the transform to each row of :math:`\mathbf{G}` and forming a new matrix consisting of rows of transformed sensitivity is equivalent to the following operation:

.. math::
     \widetilde{\mathbf{G}}=\mathbf{G}\mathcal{W}^T,
     :label: senswvt

where :math:`\widetilde{\mathbf{G}}` is the transformed matrix. The thresholding is applied to individual rows of :math:`\mathbf{G}` by the following rule to form the sparse representation :math:`\widetilde{\mathbf{G}}^S`,

.. math::
     \widetilde{g}_{ij}^{s}=\begin{cases} \widetilde{g}_{ij} & \mbox{if } \left|\widetilde{g}_{ij}\right| \geq \delta _i \\
     0 & \mbox{if } \left|\widetilde{g}_{ij}\right| < \delta _i
     \end{cases}, ~~ i=1,\ldots,N,
     :label: elemg


where :math:`\delta _i` is the threshold level, and :math:`\widetilde{g}_{ij}` and :math:`\widetilde{g}_{ij}^{s}` are the elements of :math:`\widetilde{\mathbf{G}}` and :math:`\widetilde{\mathbf{G}}^S`, respectively. The threshold level :math:`\delta _i` are determined according to the allowable error of the reconstructed sensitivity, which is measured by the ratio of norm of the error in each row to the norm of that row, :math:`r_i(\delta_i)`. It can be evaluated directly in the wavelet domain by the following expression:

.. math:: 
    r_i(\delta_i)=\sqrt{\frac{\underset{\left | {\widetilde{g}_{ij}} \right |<\delta_i}\sum{\widetilde{g}_{ij}}^2}{\underset{j}\sum{\widetilde{g}_{ij}^2}}}, ~~i=1,\ldots,N,
    :label: rhoi


Here the numerator is the norm of the discarded coefficients and the denominator is the norm of all coefficients. The threshold level :math:`\delta_{i_o}` is calculated on a representative row, :math:`i_o`. This threshold is then used to define a relative threshold :math:`\epsilon =\delta_{i_{o}}/ \underset{j}{\max}\left | {\widetilde{g}_{ij}} \right |`. The absolute threshold level for each row is obtained by

.. math::
     \delta_i = \epsilon \underset{j}{\max}\left | {\widetilde{g}_{ij}} \right|, ~~i=1,\ldots,N.
     :label: deltai

The program that implements this compression procedure is GZSENES. The user is asked to specify the relative error :math:`r^*` and the program will determine the relative threshold level :math:`\delta_i`. Usually a value of a few percent is appropriate for :math:`r^*`. When both surface and borehole data are present, two different relative threshold levels are calculated by choosing a representative row for surface data and another for borehole data. For experienced users and ones that are re-inverting the data, the program also allows the direct input of the relative threshold level.

