.. role:: cpp(code)
  :language: c++

.. _sec:Verification:

Verification
============

The following test are used to verify the correct behavior of REMORA. Every problem below has a script ``test.sh``

.. _advection:

Advection
---------

The ``Advection`` problem tests scalar advection on a doubly-periodic domain with
flat bathymetry.


.. _channeltest:

Channel Test
------------

The reentrant channel test of [paper] is reproduced. This problem tests the development of turbulence and GLS mixing scheme. The [quantity] agrees with that from ROMS.


.. _doublegyre:

Double Gyre
-----------

This reproduces the classic wind-driven double gyre problem, similar to the ROMS test problem by the same name.


.. _doublyperiodic:

Doubly Periodic
---------------

The basic version of this test simulates a flow with a depth-dependent horizontal velocity and temperature profile in a doubly-periodic domain with flat bathymetry. When non-flat bathymetry is used, the depth profile is the same as in the ROMS (and REMORA) Upwelling problem.

.. _idealminigrid:

Ideal Mini Grid
---------------

This small idealized grid is used to test netCDF-provided initial and boundary conditions. The ocean is initialized with zero velocity and a constant temperature and salinity. Time-varying boundary conditions are then applied for velocity, temperature, or salinity (provided by netCDF file). The default is to used a clamped boundary condition for all quantities, but options for Chapman-Flather and radiation conditions are available. This test also verifies correct behavior with land-sea masking when using the ``_masked`` grid file.

The netCDF files needed to run these tests can be found in the `remora-data <https://github.com/seahorce-scidac/remora-data`_ repository under the ``IdealMiniGrid`` directory.

.. _particlesseamount:

Particles Over Seamount
-----------------------

This problem tests advection of tracer particles on a flat domain.

.. _seamount-desc:

Seamount
--------

The `Seamount <https://www.myroms.org/wiki/SEAMOUNT_CASE>`_ problem involves an (analytically) stably stratified fluid at rest over a seamount. In the absence of numerical errors, the fluid will remain at rest. However, this may not occur due to numerical errors in the calculation of the horizontal pressure gradient when the vertical coordinates are misaligned with the geopotential surfaces, as is the case in problems with spatially-varying bathymetry in ROMS/REMORA.


.. _upwelling-desc:

Upwelling
---------

The `Upwelling <https://www.myroms.org/wiki/UPWELLING_CASE>`_ demonstrates wind-driven upwelling over a perioidc channel. It closely matches the test problem by the same name in ROMS.
