
 .. _RegressionTests:

Regression Tests
================

There are currently 27 tests which are run as part of every PR on a range of architectures: 20 accuracy
tests that compare a plotfile against a stored reference, two that assert values against an external
reference instead, and five that assert a misconfigured input fails for the stated reason.
The CI tests use cmake and are based on the version
of AMReX in the REMORA submodule. This suite can be run following the
instructions in :ref:`Testing<Testing>`, which also describes the different kinds of assertion and what
each is for. All of them run under MPI on two ranks, and are checked to be rank-invariant.

In addition there is a suite of more extensive nightly tests that use GNUMake and use the current
development branch of AMReX.

Results from the nightly CPU tests can be found here: `CPU tests`_

Results from the nightly GPU tests can be found here: `GPU tests`_

.. _`CPU tests`: https://ccse.lbl.gov/pub/RegressionTesting1/REMORA

.. _`GPU tests`: https://ccse.lbl.gov/pub/GpuRegressionTesting/REMORA

Continuous Integration (CI) Tests
---------------------------------

The following problems are currently tested in the CI. More details about the problems underlying these tests are given in :ref:`sec:verification`.

+-------------------------+----------+-----------+----------+---------------------------------+
| Test                    | nx ny nz | xbc       | ybc      | Other                           |
+=========================+==========+===========+==========+=================================+
| Advection               | 81 81 16 | Periodic  | Periodic |                                 |
+-------------------------+----------+-----------+----------+---------------------------------+
| Advection_ML            | 80 80 16 | Periodic  | Periodic | multilevel                      |
+-------------------------+----------+-----------+----------+---------------------------------+
| BoundaryLayer           | 39 4 30  | Radiation | Periodic | Coriolis                        |
|                         |          |           |          |                                 |
|                         |          | / outflow |          | GLS mixing scheme               |
|                         |          |           |          |                                 |
|                         |          |           |          | non-flat bathymetry             |
+-------------------------+----------+-----------+----------+---------------------------------+
| Channel_Test            | 20 60 50 | Periodic  | SlipWall | Coriolis                        |
|                         |          |           |          |                                 |
|                         |          |           |          | GLS mixing scheme               |
|                         |          |           |          |                                 |
|                         |          |           |          | non-flat bathymetry             |
|                         |          |           |          |                                 |
|                         |          |           |          | quadratic bottom stress         |
|                         |          |           |          |                                 |
|                         |          |           |          | bulk fluxes                     |
|                         |          |           |          |                                 |
|                         |          |           |          | cloud cover                     |
|                         |          |           |          |                                 |
|                         |          |           |          | evaporation/precipitation with  |
|                         |          |           |          |                                 |
|                         |          |           |          | sea surface height correction   |
+-------------------------+----------+-----------+----------+---------------------------------+
| DogboneAnalytic         | 42 15 16 | SlipWall  | SlipWall | quadratic bottom stress         |
|                         |          |           |          |                                 |
|                         |          |           |          | land-sea masking                |
+-------------------------+----------+-----------+----------+---------------------------------+
| DogboneAnalytic_MLquad  | 42 15 16 | SlipWall  | SlipWall | quadratic bottom stress         |
|                         |          |           |          |                                 |
|                         |          |           |          | land-sea masking                |
|                         |          |           |          |                                 |
|                         |          |           |          | static multilevel               |
+-------------------------+----------+-----------+----------+---------------------------------+
| DogboneAnalytic_MLvel   | 42 15 16 | SlipWall  | SlipWall | quadratic bottom stress         |
|                         |          |           |          |                                 |
|                         |          |           |          | land-sea masking                |
|                         |          |           |          |                                 |
|                         |          |           |          | dynamic multilevel              |
+-------------------------+----------+-----------+----------+---------------------------------+
| DoubleGyre              | 54 108 4 | SlipWall  | SlipWall | Coriolis                        |
+-------------------------+----------+-----------+----------+---------------------------------+
| DoublyPeriodic          | 41 80 16 | Periodic  | Periodic | Coriolis                        |
+-------------------------+----------+-----------+----------+---------------------------------+
| DoublyPeriodic_bathy    | 41 80 16 | Periodic  | Periodic | Coriolis                        |
|                         |          |           |          |                                 |
|                         |          |           |          | non-flat bathymetry             |
+-------------------------+----------+-----------+----------+---------------------------------+
| Seamount                | 49 48 13 | Periodic  | Periodic | Coriolis                        |
+-------------------------+----------+-----------+----------+---------------------------------+
| Upwelling               | 41 80 16 | Periodic  | SlipWall | Coriolis                        |
|                         |          |           |          |                                 |
|                         |          |           |          | non-flat bathymetry             |
+-------------------------+----------+-----------+----------+---------------------------------+
| Upwelling_GLS           | 41 80 16 | Periodic  | SlipWall | Coriolis                        |
|                         |          |           |          |                                 |
|                         |          |           |          | non-flat bathymetry             |
|                         |          |           |          |                                 |
|                         |          |           |          | GLS mixing scheme               |
+-------------------------+----------+-----------+----------+---------------------------------+
| Upwelling_NLEOS         | 41 80 16 | Periodic  | SlipWall | Coriolis                        |
|                         |          |           |          |                                 |
|                         |          |           |          | non-flat bathymetry             |
|                         |          |           |          |                                 |
|                         |          |           |          | nonlinear equation of state     |
+-------------------------+----------+-----------+----------+---------------------------------+
| Upwelling_logdrag       | 41 80 16 | Periodic  | SlipWall | Coriolis                        |
|                         |          |           |          |                                 |
|                         |          |           |          | non-flat bathymetry             |
|                         |          |           |          |                                 |
|                         |          |           |          | logarithmic bottom stress       |
+-------------------------+----------+-----------+----------+---------------------------------+
| Upwelling_qdrag         | 41 80 16 | Periodic  | SlipWall | Coriolis                        |
|                         |          |           |          |                                 |
|                         |          |           |          | non-flat bathymetry             |
|                         |          |           |          |                                 |
|                         |          |           |          | quadratic bottom stress         |
+-------------------------+----------+-----------+----------+---------------------------------+

.. _sec:ci-hires-tests:

High-resolution initialization tests
------------------------------------

These cover ``remora.hires_grid_level``: bathymetry evaluated on a refined level and averaged down to
level 0 (see :ref:`Inputs<sec:Inputs>`). They are analytic, so they need no NetCDF build and no external
data. ``remora.hires_init_level`` is NetCDF-only and is covered by the developer lanes in
``Exec/GulfRefinementTest`` instead, described in that directory's ``README.rst``.

The lanes come in three flavours, because "the run completed" is not evidence that an average-down
happened, let alone that it was right.

*Transparency.* ``ChannelTest`` sets ``h = 50`` and ``DogboneAnalytic`` sets ``h = 10``, both with
``setVal``, covering every grow cell. The mean over any refined block is then the same constant, so the
hires path must reproduce the corresponding non-hires baseline exactly. These lanes compare against the
existing ``Channel_Test`` and ``DogboneAnalytic_MLvel`` reference files rather than storing new ones.
They deliberately *cannot* detect a hires flag that is being ignored -- that is the next flavour's job.

*Discrimination.* ``Seamount``'s bathymetry is a Gaussian in physical position, so the mean of the
refined samples in a coarse cell is not the coarse midpoint sample. ``Seamount_hires`` must therefore
both match its own reference and **differ** from the plain ``Seamount`` one; the difference is about
16 m at the summit, some 12 orders of magnitude above the comparison tolerance. Without that second
clause a feature that silently stops taking effect still matches its own snapshot -- which is exactly
what a misspelled parameter name did in ``Exec/GulfRefinementTest`` before this was tested.

*Misconfiguration.* Five lanes assert that a bad combination fails with the message that names it, rather
than segfaulting, writing out of bounds, or running on data it quietly ignored.

+---------------------------------+----------+---------------------------------------------------------+
| Test                            | nx ny nz | What it asserts                                         |
+=================================+==========+=========================================================+
| Channel_Test_hires              | 20 60 50 | constant bathymetry averages down exactly (ratio 3),    |
|                                 |          |                                                         |
|                                 |          | compared against the Channel_Test reference             |
+---------------------------------+----------+---------------------------------------------------------+
| DogboneAnalytic_MLhires         | 42 15 16 | the same, with a live refined level and regridding,     |
|                                 |          |                                                         |
|                                 |          | so it also covers the regrid paths and masking          |
+---------------------------------+----------+---------------------------------------------------------+
| Seamount_hires                  | 49 48 13 | a Gaussian bathymetry averages down to something        |
|                                 |          |                                                         |
|                                 |          | different from the direct level-0 evaluation (ratio 3)  |
+---------------------------------+----------+---------------------------------------------------------+
| Seamount_hires_r4               | 49 48 13 | the same at ratio 4, wider than the standard grow       |
|                                 |          |                                                         |
|                                 |          | ring, and differing from the ratio-3 result             |
+---------------------------------+----------+---------------------------------------------------------+
| Upwelling_Fennel_hires_init     | 41 80 16 | biology coexists with an averaged-down bathymetry:      |
|                                 |          |                                                         |
|                                 |          | the six constant Fennel tracers and the zero dye are    |
|                                 |          |                                                         |
|                                 |          | asserted by value, catching a component-offset error    |
+---------------------------------+----------+---------------------------------------------------------+
| Seamount_hires_init_abort       | 49 48 13 | high-resolution *initialization* with analytic initial  |
|                                 |          |                                                         |
|                                 |          | conditions is rejected (it is NetCDF-only)              |
+---------------------------------+----------+---------------------------------------------------------+
| Seamount_hires_grid_max_abort   | 49 48 13 | a hires level above ``amr.max_level`` is rejected       |
+---------------------------------+----------+---------------------------------------------------------+
| Seamount_hires_init_max_abort   | 49 48 13 | the same for the initialization level                   |
+---------------------------------+----------+---------------------------------------------------------+
| Seamount_hires_grid_zero_abort  | 49 48 13 | a hires level of 0 is rejected rather than              |
|                                 |          |                                                         |
|                                 |          | dereferencing an array only allocated above level 0     |
+---------------------------------+----------+---------------------------------------------------------+
| Seamount_hires_flat_abort       | 49 48 13 | ``remora.flat_bathymetry`` combined with a hires level  |
|                                 |          |                                                         |
|                                 |          | is rejected instead of silently winning at level 0      |
+---------------------------------+----------+---------------------------------------------------------+

``Upwelling_Fennel_hires_init`` runs to zero steps: it writes the initial plotfile and exits, which takes
about a second and keeps the time stepper out of an assertion about initialization. The six Fennel
constants come from the analytic profile in ``Source/Prob/REMORA_InitAnalyticBiology_BioToy.H``, and a
constant is unchanged by averaging, so the expected level-0 values are known outright rather than
snapshotted. Asserting the dye is zero alongside them is what catches a regression in the
``Bio_comp = Tracer_comp + remora.nscalar`` offset, which would shift a biology tracer into the dye slot.

Nightly Regression Tests on CPU
-------------------------------
And the following are currently tested nighly on CPU.

Based on :ref:`Advection<advection>`:

+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Test                                   | nx ny nz     | xbc              | ybc               | Other                            |
+========================================+==============+==================+===================+==================================+
| Advection-1grid-xy                     | 81 81 16     | Periodic         | Periodic          |                                  |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Advection-1grid-xy-ML                  | 80 80 16     | Periodic         | Periodic          | multilevel                       |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Advection-OMP-xy                       | 81 81 16     | Periodic         | Periodic          | MPI + OpenMP                     |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Advection-OMP-xy-ML                    | 80 80 16     | Periodic         | Periodic          | MPI + OpenMP                     |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | multilevel                       |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Advection-xy                           | 81 81 16     | Periodic         | Periodic          | MPI                              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Advection-xy-restart                   | 81 81 16     | Periodic         | Periodic          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | restart                          |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Advection-xy-restart-ML                | 80 80 16     | Periodic         | Periodic          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | restart                          |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | multilevel                       |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Advection-U3-xy                        | 81 81 16     | Periodic         | Periodic          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | advection: upstream 3rd order    |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+

Based on :ref:`Channel Test<channeltest>`, which always includes Coriois, GLS mixing scheme, and non-flat bathymetry:

+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Test                                   | nx ny nz     | xbc              | ybc               | Other                            |
+========================================+==============+==================+===================+==================================+
| ChannelTest                            | 20 60 50     | Periodic         | SlipWall          | MPI                              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| ChannelTest-OMP                        | 20 60 50     | Periodic         | SlipWall          | MPI + OpenMP                     |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| ChannelTest-OMP-1grid-xy               | 20 60 50     | Periodic         | SlipWall          | OpenMP                           |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| ChannelTest-xy-restart                 | 20 60 50     | Periodic         | SlipWall          | MPI + OpenMP                     |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | restart                          |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| ChannelTestOrlanski                    | 20 60 50     | Radiation        | Radiation         | MPI                              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| ChannelTestOrlanski-OMP                | 20 60 50     | Radiation        | Radiation         | MPI + OpenMP                     |
|                                        |              |                  |                   |                                  |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| ChannelTestOrlanski-OMP-1grid-xy       | 20 60 50     | Radiation        | Radiation         | OpenMP                           |
|                                        |              |                  |                   |                                  |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| ChannelTestOrlanski-xy-restart         | 20 60 50     | Radiation        | Radiation         | MPI + OpenMP                     |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | restart                          |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+

Based on :ref:`Dogbone<dogbone>`, which always includes PnetCDF and land-sea masking:

+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| Test                                   | nx ny nz     | xbc              | ybc               | Other                                |
+========================================+==============+==================+===================+======================================+
| Dogbone                                | 42 15 10     | Slipwall         | SlipWall          | MPI                                  |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| Dogbone-1grid                          | 42 15 10     | Slipwall         | SlipWall          |                                      |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| Dogbone-OMP                            | 42 15 10     | Slipwall         | SlipWall          | MPI + OpenMP                         |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| Dogbone-OMP-1grid                      | 42 15 10     | Slipwall         | SlipWall          | OpenMP                               |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| Dogbone-restart                        | 42 15 10     | Slipwall         | Slipwall          | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | restart                              |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+


Based on :ref:`Double Gyre<doublegyre>`, which always includes Coriolis:

+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Test                                   | nx ny nz     | xbc              | ybc               | Other                            |
+========================================+==============+==================+===================+==================================+
| DoubleGyre                             | 54 108 4     | SlipWall         | SlipWall          | MPI                              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DoubleGyre-OMP                         | 54 108 4     | SlipWall         | SlipWall          | MPI + OpenMP                     |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DoubleGyre-OMP-1grid-xy                | 54 108 4     | SlipWall         | SlipWall          | OpenMP                           |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DoubleGyre-OMP-xy-restart              | 54 108 4     | SlipWall         | SlipWall          | MPI + OpenMP                     |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | restart                          |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+

Based on :ref:`Doubly Periodic<doublyperiodic>`, which always includes Coriolis:

+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Test                                   | nx ny nz     | xbc              | ybc               | Other                            |
+========================================+==============+==================+===================+==================================+
| DoublyPeriodic-1grid-xy                | 41 80 16     | Periodic         | Periodic          |                                  |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DoublyPeriodic-1grid-xy-bathy          | 41 80 16     | Periodic         | Periodic          | non-flat bathyemtry              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DoublyPeriodic-NETCDF-build            | 41 80 16     | N/A              | N/A               | Build w/PnetCDF                  |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DoublyPeriodic-OMP-1grid-xy            | 41 80 16     | Periodic         | Periodic          | OpenMP                           |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DoublyPeriodic-OMP-xy                  | 41 80 16     | Periodic         | Periodic          | MPI + OpenMP                     |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DoublyPeriodic-OMP-xy-bathy            | 41 80 16     | Periodic         | Periodic          | MPI + OpenMP                     |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | non-flat bathymetry              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DoublyPeriodic-xy                      | 41 80 16     | Periodic         | Periodic          | MPI                              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DoublyPeriodic-xy-bathy                | 41 80 16     | Periodic         | Periodic          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | non-flat bathymetry              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DoublyPeriodic-xy-restart              | 41 80 16     | Periodic         | Periodic          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | restart                          |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DoublyPeriodic64-OMP-xy                | 328 320 64   | Periodic         | Periodic          | MPI + OpenMP,                    |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | large problem                    |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DoublyPeriodic64-OMP-xy-bathy          | 328 320 64   | Periodic         | Periodic          | MPI + OpenMP                     |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | large problem                    |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | non-flat bathymetry              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DoublyPeriodic64-xy                    | 328 320 64   | Periodic         | Periodic          | MPI, large problem               |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DoublyPeriodicC4-xy                    | 41 80 16     | Periodic         | Periodic          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | adv.: centered 4th order         |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+

Based on :ref:`Ideal Mini Grid<idealminigrid>`, which always includes Coriolis and PnetCDF. Replace braces in test name with ``IdealMiniGrid``. C-F stands for Chapman-Flather:

+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| Test                                   | nx ny nz     | xbc              | ybc               | Other                                |
+========================================+==============+==================+===================+======================================+
| {}                                     | 10 16 20     | Clamped          | Clamped           | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary salt at boundary                |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}-1grid                               | 10 16 20     | Clamped          | Clamped           | Vary salt at boundary                |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}-CF-Uvel-OMP                         | 10 16 20     | C-F              | C-F               | MPI + OpenMP                         |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary velocity at boundary            |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}-CFO-Salt-OMP                        | 10 16 20     | C-F              | C-F               | MPI + OpenMP                         |
|                                        |              |                  |                   |                                      |
|                                        |              | Radiation        | Radiation         | Vary salt at boundary                |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}-CFO-Temp-OMP                        | 10 16 20     | C-F              | C-F               | MPI + OpenMP                         |
|                                        |              |                  |                   |                                      |
|                                        |              | Radiation        | Radiation         | Vary temperature at boundary         |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}-CFO-Uvel-OMP                        | 10 16 20     | C-F              | C-F               | MPI + OpenMP                         |
|                                        |              |                  |                   |                                      |
|                                        |              | Radiation        | Radiation         | Vary velocity at boundary            |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}-EWWall-OMP                          | 10 16 20     | SlipWall         | Clamped           | MPI + OpenMP                         |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary salt at boundary                |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}-EWWall-restart                      | 10 16 20     | SlipWall         | Clamped           | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary salt at boundary                |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | restart                              |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}-NSWall-OMP                          | 10 16 20     | Clamped          | SlipWall          | MPI + OpenMP                         |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary salt at boundary                |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}-NSWall-restart                      | 10 16 20     | Clamped          | Slipwall          | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary salt at boundary                |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | restart                              |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}-OMP                                 | 10 16 20     | Clamped          | Clamped           | MPI + OpenMP                         |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary salt at boundary                |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}-OMP-1grid                           | 10 16 20     | Clamped          | Clamped           | OpenMP                               |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary salt at boundary                |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}-Temp                                | 10 16 20     | Clamped          | Clamped           | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary temperature at boundary         |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}-Uvel                                | 10 16 20     | Clamped          | Clamped           | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary velocity at boundary            |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}-Uvel-1grid                          | 10 16 20     | Clamped          | Clamped           | Vary velocity at boundary            |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}-Uvel-EWWall-OMP                     | 10 16 20     | SlipWall         | Clamped           | MPI + OpenMP                         |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary velocity at boundary            |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}-Uvel-NSWall-OMP                     | 10 16 20     | Clamped          | SlipWall          | MPI + OpenMP                         |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary velocity at boundary            |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}-Uvel-OMP                            | 10 16 20     | Clamped          | Clamped           | MPI + OpenMP                         |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary velocity at boundary            |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}-Uvel-OMP-1grid                      | 10 16 20     | Clamped          | Clamped           | OpenMP                               |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary velocity at boundary            |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}-restart                             | 10 16 20     | Clamped          | Clamped           | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary salt at boundary                |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | restart                              |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Mask                                 | 10 16 20     | Clamped          | Clamped           | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary salt at boundary                |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | land-sea masking                     |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Mask-1grid                           | 10 16 20     | Clamped          | Clamped           | Coriolis                             |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | land-sea masking                     |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Mask-CF-Uvel-OMP                     | 10 16 20     | C-F              | C-F               | MPI + OpenMP                         |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary velocity at boundary            |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | land-sea masking                     |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Mask-CFO-Salt-OMP                    | 10 16 20     | C-F              | C-F               | MPI + OpenMP                         |
|                                        |              |                  |                   |                                      |
|                                        |              | Radiation        | Radiaion          | Vary salt at boundary                |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | land-sea masking                     |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Mask-CFO-Temp-OMP                    | 10 16 20     | C-F              | C-F               | MPI + OpenMP                         |
|                                        |              |                  |                   |                                      |
|                                        |              | Radiation        | Radiaion          | Vary temperature at boundary         |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | land-sea masking                     |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Mask-CFO-Uvel-OMP                    | 10 16 20     | C-F              | C-F               | MPI + OpenMP                         |
|                                        |              |                  |                   |                                      |
|                                        |              | Radiation        | Radiaion          | Vary velocity at boundary            |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | land-sea masking                     |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Mask-EWWall-OMP                      | 10 16 20     | SlipWall         | Clamped           | MPI + OpenMP                         |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary salt at boundary                |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | land-sea masking                     |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Mask-OMP                             | 10 16 20     | Clamped          | Clamped           | MPI + OpenMP                         |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary salt at boundary                |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | land-sea masking                     |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Mask-OMP-1grid                       | 10 16 20     | Clamped          | Clamped           | OpenMP                               |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary salt at boundary                |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | land-sea masking                     |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Mask-Temp                            | 10 16 20     | Clamped          | Clamped           | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary temperature at boundary         |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | land-sea masking                     |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Mask-Uvel                            | 10 16 20     | Clamped          | Clamped           | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary velocity at boundary            |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | land-sea masking                     |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Mask-Uvel-1grid                      | 10 16 20     | Clamped          | Clamped           | Vary velocity at boundary            |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | land-sea masking                     |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Mask-Uvel-EWWall-OMP                 | 10 16 20     | SlipWall         | Clamped           | MPI + OpenMP                         |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary velocity at boundary            |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | land-sea masking                     |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Mask-Uvel-NSWall-OMP                 | 10 16 20     | Clamped          | SlipWall          | MPI + OpenMP                         |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary velocity at boundary            |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | land-sea masking                     |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Mask-Uvel-OMP                        | 10 16 20     | Clamped          | Clamped           | MPI + OpenMP                         |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary velocity at boundary            |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | land-sea masking                     |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Mask-Uvel-OMP-1grid                  | 10 16 20     | Clamped          | Clamped           | OpenMP                               |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary velocity at boundary            |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | land-sea masking                     |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Mask-restart                         | 10 16 20     | Clamped          | Clamped           | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary salt at boundary                |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | land-sea masking                     |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | restart                              |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Wind                                 | 10 16 20     | Clamped          | Clamped           | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary salt at boundary                |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Surface wind and bulk fluxes         |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Wind-1grid                           | 10 16 20     | Clamped          | Clamped           | Vary salt at boundary                |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Surface wind and bulk fluxes         |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Wind-OMP                             | 10 16 20     | Clamped          | Clamped           | MPI + OpenMP                         |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary salt at boundary                |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Surface wind and bulk fluxes         |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Wind-OMP-1grid                       | 10 16 20     | Clamped          | Clamped           | OpenMP                               |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary salt at boundary                |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Surface wind and bulk fluxes         |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Wind-restart                         | 10 16 20     | Clamped          | Clamped           | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary salt at boundary                |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | restart                              |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Surface wind and bulk fluxes         |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Clim                                 | 10 16 20     | C-F              | C-F               | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              | Radiation        | Radiation         | Vary salt, temp, v at boundary       |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Climatology nudging                  |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Clim-OMP                             | 10 16 20     | C-F              | C-F               | MPI + OpenMP                         |
|                                        |              |                  |                   |                                      |
|                                        |              | Radiation        | Radiation         | Vary salt, temp, v at boundary       |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Climatology nudging                  |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Clim-OMP-1grid                       | 10 16 20     | C-F              | C-F               | OpenMP                               |
|                                        |              |                  |                   |                                      |
|                                        |              | Radiation        | Radiation         | Vary salt, temp, v at boundary       |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Climatology nudging                  |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Clim-restart                         | 10 16 20     | C-F              | C-F               | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              | Radiation        | Radiation         | Vary salt, temp, v at boundary       |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | restart                              |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Climatology nudging                  |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+

Based on :ref:`Ideal River Grid<idealrivgrid>`, which always includes Coriolis, PnetCDF, and rivers:

+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| Test                                   | nx ny nz     | xbc              | ybc               | Other                                |
+========================================+==============+==================+===================+======================================+
| IdealRivGrid                           | 10 16 20     | Clamped          | Clamped           | MPI                                  |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| IdealRivGrid-1grid                     | 10 16 20     | Clamped          | Clamped           |                                      |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| IdealRivGrid-OMP                       | 10 16 20     | Clamped          | Clamped           | MPI + OpenMP                         |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| IdealRivGrid-OMP-1grid                 | 10 16 20     | Clamped          | Clamped           | OpenMP                               |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| IdealRivGrid-restart                   | 10 16 20     | Clamped          | Clamped           | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | restart                              |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+

Based on :ref:`Particles Over Seamount<particlesseamount>`, which always include MPI, Coriolis, and tracer particles:

+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Test                                   | nx ny nz     | xbc              | ybc               | Other                            |
+========================================+==============+==================+===================+==================================+
| ParticlesOverSeamount                  | 41 80 16     | Periodic         | Periodic          |                                  |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| ParticlesOverSeamount-restart          | 41 80 16     | Periodic         | Periodic          | restart                          |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+

Based on :ref:`Seamount<seamount-desc>`, which always includes Coriolis and non-flat bathymetry:

+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Test                                   | nx ny nz     | xbc              | ybc               | Other                            |
+========================================+==============+==================+===================+==================================+
| Seamount-1grid-xy                      | 49 48 13     | Periodic         | Periodic          |                                  |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Seamount-OMP-xy                        | 49 48 13     | Periodic         | Periodic          | MPI + OpenMP                     |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Seamount-xy                            | 49 48 13     | Periodic         | Periodic          | MPI                              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Seamount-xy-restart                    | 49 48 13     | Periodic         | Periodic          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | restart                          |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Seamount64-OMP-xy                      | 320 320 64   | Periodic         | Periodic          | MPI + OpenMP, large problem      |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+

Based on :ref:`Upwelling<upwelling-desc>`, which always includes Coriolis and non-flat bathymetry:

+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Test                                   | nx ny nz     | xbc              | ybc               | Other                            |
+========================================+==============+==================+===================+==================================+
| Upwelling                              | 41 80 16     | Periodic         | SlipWall          | MPI                              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling-1grid                        | 41 80 16     | Periodic         | SlipWall          | Coriolis                         |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling-OMP                          | 41 80 16     | Periodic         | SlipWall          | MPI + OpenMP                     |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling-logDrag-OMP                  | 41 80 16     | Periodic         | SlipWall          | MPI + OpenMP                     |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | logarithmic bottom stress        |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling-nonlinEOS-OMP                | 41 80 16     | Periodic         | SlipWall          | MPI + OpenMP                     |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | nonlinear equation of state      |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling-quadDrag-OMP                 | 41 80 16     | Periodic         | SlipWall          | MPI + OpenMP                     |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | quadratic bottom stress          |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling-OMP-1grid                    | 41 80 16     | Periodic         | SlipWall          | OpenMP                           |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling-restart                      | 41 80 16     | Periodic         | SlipWall          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | restart                          |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling-x                            | 41 80 16     | SlipWall         | Periodic          | MPI                              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling-x-1grid                      | 41 80 16     | SlipWall         | Periodic          |                                  |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling-x-OMP                        | 41 80 16     | SlipWall         | Periodic          | MPI + OpenMP                     |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling64-OMP                        | 328 320 64   | SlipWall         | Periodic          | MPI + OpenMP, large problem      |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling64-OMP                        | 328 320 64   | SlipWall         | Periodic          | MPI + OpenMP, large problem      |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| UpwellingC4                            | 41 80 16     | Periodic         | SlipWall          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | advection: centered 4th order    |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling_GLS                          | 41 80 16     | Periodic         | SlipWall          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | GLS mixing scheme                |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling_GLS-restart                  | 41 80 16     | Periodic         | SlipWall          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | GLS mixing scheme                |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | restart                          |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling_GLS_Canuto_A                 | 41 80 16     | Periodic         | SlipWall          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | GLS mixing scheme                |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | Canuto A stability               |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling_GLS_Canuto_B                 | 41 80 16     | Periodic         | SlipWall          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | GLS mixing scheme                |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | Canuto B stability               |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+

Nightly Regression Tests on GPU
-------------------------------

And the following are currently tested nighly on GPU. All are compiled and run with CUDA.

Based on :ref:`Advection<advection>`:

+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Test                                   | nx ny nz     | xbc              | ybc               | Other                            |
+========================================+==============+==================+===================+==================================+
| Advection-1grid-xy                     | 81 81 16     | Periodic         | Periodic          |                                  |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Advection-1grid-xy-ML                  | 80 80 16     | Periodic         | Periodic          | multilevel                       |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Advection-xy                           | 81 81 16     | Periodic         | Periodic          | MPI                              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Advection-xy-ML                        | 80 80 16     | Periodic         | Periodic          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | multilevel                       |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Advection-xy-restart                   | 81 81 16     | Periodic         | Periodic          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | restart                          |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Advection64-xy                         | 328 320 64   | Periodic         | Periodic          | MPI, large problem               |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | restart                          |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Advection-U3-xy                        | 81 81 16     | Periodic         | Periodic          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | advection: upstream 3rd order    |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+

Based on :ref:`Channel Test<channeltest>`, which always includes Coriolis, GLS mixing scheme, and non-flat bathymetry:

+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Test                                   | nx ny nz     | xbc              | ybc               | Other                            |
+========================================+==============+==================+===================+==================================+
| ChannelTest-1grid-xy                   | 20 60 50     | Periodic         | SlipWall          |                                  |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| ChannelTest-xy                         | 20 60 50     | Periodic         | SlipWall          | MPI                              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| ChannelTest-xy-restart                 | 20 60 50     | Periodic         | SlipWall          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | restart                          |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+

Based on :ref:`Dogbone<dogbone>`, which always includes PnetCDF and land-sea masking:

+----------------------------------------+--------------+------------------+-------------------+-------------+
| Test                                   | nx ny nz     | xbc              | ybc               | Other       |
+========================================+==============+==================+===================+=============+
| Dogbone                                | 42 15 10     | Slipwall         | Slipwall          | MPI         |
+----------------------------------------+--------------+------------------+-------------------+-------------+
| Dogbone-1grid                          | 42 15 10     | Slipwall         | Slipwall          |             |
+----------------------------------------+--------------+------------------+-------------------+-------------+
| Dogbone-restart                        | 42 15 10     | Slipwall         | Slipwall          | MPI         |
|                                        |              |                  |                   |             |
|                                        |              |                  |                   | restart     |
+----------------------------------------+--------------+------------------+-------------------+-------------+


Based on :ref:`Double Gyre<doublegyre>`, which always includes Coriolis:

+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Test                                   | nx ny nz     | xbc              | ybc               | Other                            |
+========================================+==============+==================+===================+==================================+
| DoubleGyre-1grid-xy                    | 54 108 4     | SlipWall         | SlipWall          |                                  |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DoubleGyre-xy                          | 54 108 4     | SlipWall         | SlipWall          | MPI                              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DoubleGyre-xy-restart                  | 54 108 4     | SlipWall         | SlipWall          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | restart                          |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+

Based on :ref:`Doubly Periodic<doublyperiodic>`, which always includes Coriolis:

+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Test                                   | nx ny nz     | xbc              | ybc               | Other                            |
+========================================+==============+==================+===================+==================================+
| DoublyPeriodic-1grid-xy                | 41 80 16     | Periodic         | Periodic          |                                  |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DoublyPeriodic-xy                      | 41 80 16     | Periodic         | Periodic          | MPI                              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DoublyPeriodic-xy-bathy                | 41 80 16     | Periodic         | Periodic          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | non-flat bathymetry              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DoublyPeriodic-xy-restart              | 41 80 16     | Periodic         | Periodic          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | restart                          |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DoublyPeriodic64-xy                    | 328 320 64   | Periodic         | Periodic          | MPI, large problem               |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DoublyPeriodic64-xy-bathy              | 328 320 64   | Periodic         | Periodic          | MPI, large problem               |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | non-flat bathymetry              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DoublyPeriodicC4-xy                    | 41 80 16     | Periodic         | Periodic          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | adv.: centered 4th order         |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+

Based on :ref:`Ideal Mini Grid<idealminigrid>`, which always includes Coriolis and PnetCDF. Replace braces in test name with ``IdealMiniGrid``. C-F stands for Chapman-Flather:

+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| Test                                   | nx ny nz     | xbc              | ybc               | Other                                |
+========================================+==============+==================+===================+======================================+
| {}                                     | 10 16 20     | Clamped          | Clamped           | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary salt at boundary                |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}-1grid                               | 10 16 20     | Clamped          | Clamped           | Vary salt at boundary                |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}-CF-Uvel                             | 10 16 20     | C-F              | C-F               | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary velocity at boundary            |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}-CFO-Salt                            | 10 16 20     | C-F              | C-F               | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              | Radiation        | Radiation         | Vary salt at boundary                |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}-CFO-Temp                            | 10 16 20     | C-F              | C-F               | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              | Radiation        | Radiation         | Vary temperature at boundary         |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}-CFO-Uvel                            | 10 16 20     | C-F              | C-F               | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              | Radiation        | Radiation         | Vary velocity at boundary            |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}-EWWall                              | 10 16 20     | SlipWall         | Clamped           | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary salt at boundary                |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}-NSWall                              | 10 16 20     | Clamped          | SlipWall          | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary salt at boundary                |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}-Temp                                | 10 16 20     | Clamped          | Clamped           | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary temperature at boundary         |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}-Uvel                                | 10 16 20     | Clamped          | Clamped           | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary velocity at boundary            |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}-Uvel-1grid                          | 10 16 20     | Clamped          | Clamped           | Vary velocity at boundary            |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}-Uvel-EWWall                         | 10 16 20     | SlipWall         | Clamped           | Vary velocity at boundary            |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}-Uvel-NSWall                         | 10 16 20     | Clamped          | SlipWall          | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary velocity at boundary            |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}-restart                             | 10 16 20     | Clamped          | Clamped           | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary salt at boundary                |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | restart                              |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Mask                                 | 10 16 20     | Clamped          | Clamped           | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary salt at boundary                |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | land-sea masking                     |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Mask-1grid                           | 10 16 20     | Clamped          | Clamped           | Vary salt at boundary                |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | land-sea masking                     |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Mask-CF-Uvel                         | 10 16 20     | C-F              | C-F               | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary velocity at boundary            |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | land-sea masking                     |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Mask-CFO-Salt                        | 10 16 20     | C-F              | C-F               | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              | Radiation        | Radiaion          | Vary salt at boundary                |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | land-sea masking                     |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Mask-CFO-Temp                        | 10 16 20     | C-F              | C-F               | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              | Radiation        | Radiaion          | Vary temperature at boundary         |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | land-sea masking                     |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Mask-CFO-Uvel                        | 10 16 20     | C-F              | C-F               | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              | Radiation        | Radiaion          | Vary velocity at boundary            |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | land-sea masking                     |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Mask-EWWall                          | 10 16 20     | SlipWall         | Clamped           | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary salt at boundary                |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | land-sea masking                     |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Mask-EWWall-restart                  | 10 16 20     | SlipWall         | Clamped           | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary salt at boundary                |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | land-sea masking                     |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | restart                              |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Mask-NSWall                          | 10 16 20     | Clamped          | Slipwall          | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary salt at boundary                |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | land-sea masking                     |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Mask-NSWall-restart                  | 10 16 20     | Clamped          | Slipwall          | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary salt at boundary                |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | land-sea masking                     |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | restart                              |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Mask-Temp                            | 10 16 20     | Clamped          | Clamped           | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary temperature at boundary         |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | land-sea masking                     |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Mask-Uvel                            | 10 16 20     | Clamped          | Clamped           | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary velocity at boundary            |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | land-sea masking                     |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Mask-Uvel-1grid                      | 10 16 20     | Clamped          | Clamped           | Vary velocity at boundary            |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | land-sea masking                     |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Mask-Uvel-EWWall                     | 10 16 20     | SlipWall         | Clamped           | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary salt at boundary                |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | land-sea masking                     |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Mask-Uvel-NSWall                     | 10 16 20     | Clamped          | Periodic          | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary salt at boundary                |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | land-sea masking                     |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Mask-restart                         | 10 16 20     | Clamped          | Clamped           | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary salt at boundary                |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | land-sea masking                     |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | restart                              |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Wind                                 | 10 16 20     | Clamped          | Clamped           | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary salt at boundary                |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Surface wind and bulk fluxes         |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Wind-1grid                           | 10 16 20     | Clamped          | Clamped           | Vary salt at boundary                |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Surface wind and bulk fluxes         |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Wind-restart                         | 10 16 20     | Clamped          | Clamped           | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary salt at boundary                |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | restart                              |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Surface wind and bulk fluxes         |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Clim                                 | 10 16 20     | C-F              | C-F               | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              | Radiation        | Radiation         | Vary salt, temp, v at boundary       |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Climatology nudging                  |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Clim-1grid                           | 10 16 20     | C-F              | C-F               | Vary salt, temp, v at boundary       |
|                                        |              |                  |                   |                                      |
|                                        |              | Radiation        | Radiation         | Climatology nudging                  |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Clim-restart                         | 10 16 20     | C-F              | C-F               | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              | Radiation        | Radiation         | Vary salt, temp, v at boundary       |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | restart                              |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Climatology nudging                  |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+

Based on :ref:`Ideal River Grid<idealrivgrid>`, which always includes Coriolis, PnetCDF, and rivers:

+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| Test                                   | nx ny nz     | xbc              | ybc               | Other                                |
+========================================+==============+==================+===================+======================================+
| IdealRivGrid                           | 10 16 20     | Clamped          | Clamped           | MPI                                  |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| IdealRivGrid-1grid                     | 10 16 20     | Clamped          | Clamped           |                                      |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| IdealRivGrid-restart                   | 10 16 20     | Clamped          | Clamped           | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | restart                              |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+

Based on :ref:`Particles Over Seamount<particlesseamount>`, which always includes MPI, Coriolis, and tracer particles:

+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Test                                   | nx ny nz     | xbc              | ybc               | Other                            |
+========================================+==============+==================+===================+==================================+
| ParticlesOverSeamount                  | 41 80 16     | Periodic         | Periodic          |                                  |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| ParticlesOverSeamount-restart          | 41 80 16     | Periodic         | Periodic          | restart                          |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+

Based on :ref:`Seamount<seamount-desc>`, which always includes Coriolis and non-flat bathymetry:

+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Test                                   | nx ny nz     | xbc              | ybc               | Other                            |
+========================================+==============+==================+===================+==================================+
| Seamount-1grid-xy                      | 49 48 13     | Periodic         | Periodic          |                                  |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Seamount-xy                            | 49 48 13     | Periodic         | Periodic          | MPI                              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Seamount64-xy                          | 320 320 64   | Periodic         | Periodic          | MPI, large problem               |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+

Based on :ref:`Upwelling<upwelling-desc>`, which always includes Coriolis and non-flat bathymetry:

+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Test                                   | nx ny nz     | xbc              | ybc               | Other                            |
+========================================+==============+==================+===================+==================================+
| Upwelling                              | 41 80 16     | Periodic         | SlipWall          | MPI                              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling-1grid                        | 41 80 16     | Periodic         | SlipWall          |                                  |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling-logDrag                      | 41 80 16     | Periodic         | SlipWall          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | logarithmic bottom stress        |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling-nonlinEOS                    | 41 80 16     | Periodic         | SlipWall          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | nonlinear equation of state      |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling-quadDrag                     | 41 80 16     | Periodic         | SlipWall          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | quadratic bottom stress          |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling-x                            | 41 80 16     | SlipWall         | Periodic          | MPI                              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling-x-1grid                      | 41 80 16     | SlipWall         | Periodic          |                                  |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling64                            | 328 320 64   | SlipWall         | Periodic          | MPI, large problem               |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| UpwellingC4                            | 41 80 16     | Periodic         | SlipWall          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | advection: centered 4th order    |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling_GLS                          | 41 80 16     | Periodic         | SlipWall          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | GLS mixing scheme                |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling_GLS-restart                  | 41 80 16     | Periodic         | SlipWall          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | GLS mixing scheme                |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | restart                          |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling_GLS_Canuto_A                 | 41 80 16     | Periodic         | SlipWall          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | GLS mixing scheme                |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | Canuto A stability               |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling_GLS_Canuto_B                 | 41 80 16     | Periodic         | SlipWall          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | GLS mixing scheme                |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | Canuto B stability               |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+

