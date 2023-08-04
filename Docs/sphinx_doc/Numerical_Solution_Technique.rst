
 .. role:: cpp(code)
    :language: c++

.. _Numerical_Solution:

Numerical Solution Technique
============================
.. _Horizontal_Discretization:

Horizontal Discretization
-------------------------
In the horizontal, the ROMSX governing equations are discretized over a boundary-fitted, orthogonal curvilinear coordinates :math:`\left( \xi , \eta \right)` grid. The general formulation of the curvilinear coordinates system allows Cartesian, polar, and spherical coordinates applications. The transformation of any of these coordinates to ROMSX :math:`\left( \xi , \eta \right)` grid is specified in the metric terms (``pm``, ``pn``).

The model state variables are staggered using an ``Arakawa C-grid``. As illustrated below, the free-surface (``zeta``), density (``rho``), and active/passive tracers (``t``) are located at the center of the cell whereas the horizontal velocity (``u`` and ``v``) are located at the west/east and south/north edges of the cell, respectively. That is, the density is evaluated between points where the currents are evaluated.

Staggered Horizontal Grid
~~~~~~~~~~~~~~~~~~~~~~~~~
.. image:: figures/staggered_grid_rho_cells.png
   :width: 100%

In ROMSX all the state arrays are dimensioned the same size to facilitate parallelization. However, the computational ranges for all the state variables are:

Grid Cell
~~~~~~~~~
.. image:: figures/grid_cell.png
   :width: 50%

+---------------------------+---------------------------+-------------------------+
| Variable                  | Interior Range            | Full Range              |
+===========================+===========================+=========================+
| :math:`\rho \text{-type}` | 1:``Lm(ng)``,1:``Mm(ng)`` | 0:``L(ng)``,0:``M(ng)`` |
+---------------------------+---------------------------+-------------------------+
| :math:`\psi \text{-type}` | 2:``Lm(ng)``,2:``Mm(ng)`` | 1:``L(ng)``,1:``M(ng)`` |
+---------------------------+---------------------------+-------------------------+
| :math:`\text{u-type}`     | 2:``Lm(ng)``,1:``Mm(ng)`` | 1:``L(ng)``,0:``M(ng)`` |
+---------------------------+---------------------------+-------------------------+
| :math:`\text{v-type}`     | 1:``Lm(ng)``,2:``Mm(ng)`` | 0:``L(ng)``,1:``M(ng)`` |
+---------------------------+---------------------------+-------------------------+

.. _Vertical_Discretization:

Vertical Discretization
-----------------------
The ROMSX governing equations are discretized over variable topography using a stretched, terrain-following, vertical coordinate. As a result, each grid cell may have different level thickness (``Hz``) and volume. The model state variables are vertically staggered so that horizontal momentum (``u``, ``v``), (``rho``), and active/passive tracers (``t``) are located at the center of the grid cell. The vertical velocity (``omega``, ``w``) and vertical mixing variables (``Akt``, ``Akv``, etc) are located at the bottom and top faces of the cell. See diagram below.

Vieste-Dubrovnik Transect
~~~~~~~~~~~~~~~~~~~~~~~~~
.. image:: figures/vieste-dubrovnik.png
   :width: 100%

Staggered Vertical Grid
~~~~~~~~~~~~~~~~~~~~~~~
.. image:: figures/vertical_grid.png
   :width: 100%

The total thickness of the water column is :math:`\zeta \left( i,j\right) + h\left( i,j\right)`. The bathymetry (``h``) is usually time invariant whereas the free-surface (``zeta``) evolves in time. At input and output, the bathymetry is always a positive quantity. However, the depths ``z_r(i,j,k)`` and ``z_w(i,j,k)`` are negative for all locations below the mean sea level.

