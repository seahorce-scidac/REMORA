
 .. _RegressionTests:

Regression Tests
================

There are currently 9 accuracy tests which are run as part of every PR.
The CI tests use cmake and are based on the version
of AMReX in the REMORA submodule. This suite can be run following the
instructions in :ref:`Testing`<Testing>`.

In addition there is a suite of more extensive nightly tests that use GNUMake and use the current
development branch of AMReX.

Results from the nightly CPU tests can be found here: `CPU tests`_

Results from the nightly GPU tests can be found here: `GPU tests`_

.. _`CPU tests`: https://ccse.lbl.gov/pub/RegressionTesting1/REMORA

.. _`GPU tests`: https://ccse.lbl.gov/pub/GpuRegressionTesting/REMORA

Continuous Integration (CI) Tests
---------------------------------

The following problems are currently tested in the CI. More details about the problems underlying these tests are given in :ref:`sec:verification`.

+----------------------+----------+----------+----------+-----------------------+
| Test                 | nx ny nz | xbc      | ybc      | Other                 |
+======================+==========+==========+==========+=======================+
| Advection            | 81 81 16 | Periodic | Periodic |                       |
+----------------------+----------+----------+----------+-----------------------+
| Advection_ML         | 80 80 16 | Periodic | Periodic | multilevel            |
+----------------------+----------+----------+----------+-----------------------+
| Channel_Test         | 20 60 50 | Periodic | SlipWall | Coriolis              |
|                      |          |          |          |                       |
|                      |          |          |          | GLS mixing scheme     |
|                      |          |          |          |                       |
|                      |          |          |          | non-flat bathymetry   |
+----------------------+----------+----------+----------+-----------------------+
| DoubleGyre           | 54 108 4 | SlipWall | SlipWall | Coriolis              |
+----------------------+----------+----------+----------+-----------------------+
| DoublyPeriodic       | 41 80 16 | Periodic | Periodic | Coriolis              |
+----------------------+----------+----------+----------+-----------------------+
| DoublyPeriodic_bathy | 41 80 16 | Periodic | Periodic | Coriolis              |
|                      |          |          |          |                       |
|                      |          |          |          | non-flat bathymetry   |
+----------------------+----------+----------+----------+-----------------------+
| Seamount             | 49 48 13 | Periodic | Periodic | Coriolis              |
+----------------------+----------+----------+----------+-----------------------+
| Upwelling            | 41 80 16 | Periodic | SlipWall | Coriolis              |
|                      |          |          |          |                       |
|                      |          |          |          | non-flat bathymetry   |
+----------------------+----------+----------+----------+-----------------------+
| Upwelling_GLS        | 41 80 16 | Periodic | SlipWall | Coriolis              |
|                      |          |          |          |                       |
|                      |          |          |          | non-flat bathymetry   |
|                      |          |          |          |                       |
|                      |          |          |          | GLS mixing scheme     |
+----------------------+----------+----------+----------+-----------------------+

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
| DoublyPeriodic64-OMP-xy                | 328 320 64   | Periodic         | Periodic          | MPI + OpenMP, large problem      |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DoublyPeriodic64-OMP-xy-bathy          | 328 320 64   | Periodic         | Periodic          | MPI + OpenMP, large problem      |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | non-flat bathymetry              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DoublyPeriodic64-xy                    | 328 320 64   | Periodic         | Periodic          | MPI, large problem               |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DoublyPeriodicC4-xy                    | 41 80 16     | Periodic         | Periodic          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | advection: centered 4th order    |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+

Based on :ref:`Ideal Mini Grid<idealminigrid>`, which always includes Coriolis and PnetCDF:

+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Test                                   | nx ny nz     | xbc              | ybc               | Other                            |
+========================================+==============+==================+===================+==================================+
| IdealMiniGrid                          | 10 16 20     | Clamped          | Clamped           | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | Varying salt at boundary         |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGrid-1grid                    | 10 16 20     | Clamped          | Clamped           | Varying salt at boundary         |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGrid-CF-Uvel-OMP              | 10 16 20     | Chapman-Flather  | Chapman-Flather   | MPI + OpenMP                     |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | Varying velocity at boundary     |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGrid-CFO-Salt-OMP             | 10 16 20     | Chapman-Flather  | Chapman-Flather   | MPI + OpenMP                     |
|                                        |              |                  |                   |                                  |
|                                        |              | Radiation        | Radiation         | Varying salt at boundary         |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGrid-CFO-Temp-OMP             | 10 16 20     | Chapman-Flather  | Chapman-Flather   | MPI + OpenMP                     |
|                                        |              |                  |                   |                                  |
|                                        |              | Radiation        | Radiation         | Varying temperature at boundary  |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGrid-CFO-Uvel-OMP             | 10 16 20     | Chapman-Flather  | Chapman-Flather   | MPI + OpenMP                     |
|                                        |              |                  |                   |                                  |
|                                        |              | Radiation        | Radiation         | Varying velocity at boundary     |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGrid-EWWall-OMP               | 10 16 20     | SlipWall         | Clamped           | MPI + OpenMP                     |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | Varying salt at boundary         |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGrid-EWWall-restart           | 10 16 20     | SlipWall         | Clamped           | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | Varying salt at boundary         |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | restart                          |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGrid-NSWall-OMP               | 10 16 20     | Clamped          | SlipWall          | MPI + OpenMP                     |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | Varying salt at boundary         |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGrid-NSWall-restart           | 10 16 20     | Clamped          | Slipwall          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | Varying salt at boundary         |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | restart                          |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGrid-OMP                      | 10 16 20     | Clamped          | Clamped           | MPI + OpenMP                     |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | Varying salt at boundary         |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGrid-OMP-1grid                | 10 16 20     | Clamped          | Clamped           | OpenMP                           |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | Varying salt at boundary         |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGrid-Temp                     | 10 16 20     | Clamped          | Clamped           | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | Varying temperature at boundary  |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGrid-Uvel                     | 10 16 20     | Clamped          | Clamped           | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | Varying velocity at boundary     |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGrid-Uvel-1grid               | 10 16 20     | Clamped          | Clamped           | Varying velocity at boundary     |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGrid-Uvel-EWWall-OMP          | 10 16 20     | SlipWall         | Clamped           | MPI + OpenMP                     |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | Varying velocity at boundary     |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGrid-Uvel-NSWall-OMP          | 10 16 20     | Clamped          | SlipWall          | MPI + OpenMP                     |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | Varying velocity at boundary     |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGrid-Uvel-OMP                 | 10 16 20     | Clamped          | Clamped           | MPI + OpenMP                     |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | Varying velocity at boundary     |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGrid-Uvel-OMP-1grid           | 10 16 20     | Clamped          | Clamped           | OpenMP                           |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | Varying velocity at boundary     |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGrid-restart                  | 10 16 20     | Clamped          | Clamped           | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | Varying salt at boundary         |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | restart                          |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGridMask                      | 10 16 20     | Clamped          | Clamped           | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | Varying salt at boundary         |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | land-sea masking                 |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGridMask-1grid                | 10 16 20     | Clamped          | Clamped           | Coriolis                         |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | land-sea masking                 |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGridMask-CF-Uvel-OMP          | 10 16 20     | Chapman-Flather  | Chapman-Flather   | MPI + OpenMP                     |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | Varying velocity at boundary     |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | land-sea masking                 |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGridMask-CFO-Salt-OMP         | 10 16 20     | Chapman-Flather  | Chapman-Flather   | MPI + OpenMP                     |
|                                        |              |                  |                   |                                  |
|                                        |              | Radiation        | Radiaion          | Varying salt at boundary         |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | land-sea masking                 |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGridMask-CFO-Temp-OMP         | 10 16 20     | Chapman-Flather  | Chapman-Flather   | MPI + OpenMP                     |
|                                        |              |                  |                   |                                  |
|                                        |              | Radiation        | Radiaion          | Varying temperature at boundary  |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | land-sea masking                 |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGridMask-CFO-Uvel-OMP         | 10 16 20     | Chapman-Flather  | Chapman-Flather   | MPI + OpenMP                     |
|                                        |              |                  |                   |                                  |
|                                        |              | Radiation        | Radiaion          | Varying velocity at boundary     |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | land-sea masking                 |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGridMask-EWWall-OMP           | 10 16 20     | SlipWall         | Clamped           | MPI + OpenMP                     |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | Varying salt at boundary         |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | land-sea masking                 |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGridMask-OMP                  | 10 16 20     | Clamped          | Clamped           | MPI + OpenMP                     |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | Varying salt at boundary         |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | land-sea masking                 |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGridMask-OMP-1grid            | 10 16 20     | Clamped          | Clamped           | OpenMP                           |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | Varying salt at boundary         |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | land-sea masking                 |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGridMask-Temp                 | 10 16 20     | Clamped          | Clamped           | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | Varying temperature at boundary  |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | land-sea masking                 |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGridMask-Uvel                 | 10 16 20     | Clamped          | Clamped           | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | Varying velocity at boundary     |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | land-sea masking                 |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGridMask-Uvel-1grid           | 10 16 20     | Clamped          | Clamped           | Varying velocity at boundary     |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | land-sea masking                 |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGridMask-Uvel-EWWall-OMP      | 10 16 20     | SlipWall         | Clamped           | MPI + OpenMP                     |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | Varying velocity at boundary     |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | land-sea masking                 |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGridMask-Uvel-NSWall-OMP      | 10 16 20     | Clamped          | SlipWall          | MPI + OpenMP                     |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | Varying velocity at boundary     |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | land-sea masking                 |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGridMask-Uvel-OMP             | 10 16 20     | Clamped          | Clamped           | MPI + OpenMP                     |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | Varying velocity at boundary     |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | land-sea masking                 |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGridMask-Uvel-OMP-1grid       | 10 16 20     | Clamped          | Clamped           | OpenMP                           |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | Varying velocity at boundary     |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | land-sea masking                 |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGridMask-restart              | 10 16 20     | Clamped          | Clamped           | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | Varying salt at boundary         |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | land-sea masking                 |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | restart                          |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+

Based on :ref:`Particles Over Seamount<particlesseamount>`, which always include MPI, Coriolis, and tracer particles:

+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Test                                   | nx ny nz     | xbc              | ybc               | Other                            |
+========================================+==============+==================+===================+==================================+
| ParticlesOverSeamount                  | 41 80 16     | Periodic         | Periodic          |                                  |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| ParticlesOverSeamount-restart          | 41 80 16     | Periodic         | Periodic          | restart                          |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+

Based on :ref:`Seamount<seamount>`, which always includes Coriolis and non-flat bathymetry:

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

Based on :ref:`Upwelling<upwelling>`, which always includes Coriolis and non-flat bathymetry:

+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Test                                   | nx ny nz     | xbc              | ybc               | Other                            |
+========================================+==============+==================+===================+==================================+
| Upwelling                              | 41 80 16     | Periodic         | SlipWall          | MPI                              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling-1grid                        | 41 80 16     | Periodic         | SlipWall          | Coriolis                         |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling-OMP                          | 41 80 16     | Periodic         | SlipWall          | MPI + OpenMP                     |
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

Based on :ref:`Advection`<advection>`:

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
|                                        |              |                  |                   | advection: centered 4th order    |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+

Based on :ref:`Ideal Mini Grid<idealminigrid>`, which always includes Coriolis and PnetCDF:

+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Test                                   | nx ny nz     | xbc              | ybc               | Other                            |
+========================================+==============+==================+===================+==================================+
| IdealMiniGrid                          | 10 16 20     | Clamped          | Clamped           | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | Varying salt at boundary         |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGrid-1grid                    | 10 16 20     | Clamped          | Clamped           | Varying salt at boundary         |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGrid-CF-Uvel                  | 10 16 20     | Chapman-Flather  | Chapman-Flather   | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | Varying velocity at boundary     |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGrid-CFO-Salt                 | 10 16 20     | Chapman-Flather  | Chapman-Flather   | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              | Radiation        | Radiation         | Varying salt at boundary         |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGrid-CFO-Temp                 | 10 16 20     | Chapman-Flather  | Chapman-Flather   | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              | Radiation        | Radiation         | Varying temperature at boundary  |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGrid-CFO-Uvel                 | 10 16 20     | Chapman-Flather  | Chapman-Flather   | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              | Radiation        | Radiation         | Varying velocity at boundary     |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGrid-EWWall                   | 10 16 20     | SlipWall         | Clamped           | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | Varying salt at boundary         |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGrid-NSWall                   | 10 16 20     | Clamped          | SlipWall          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | Varying salt at boundary         |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGrid-Temp                     | 10 16 20     | Clamped          | Clamped           | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | Varying temperature at boundary  |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGrid-Uvel                     | 10 16 20     | Clamped          | Clamped           | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | Varying velocity at boundary     |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGrid-Uvel-1grid               | 10 16 20     | Clamped          | Clamped           | Varying velocity at boundary     |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGrid-Uvel-EWWall              | 10 16 20     | SlipWall         | Clamped           | Varying velocity at boundary     |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGrid-Uvel-NSWall              | 10 16 20     | Clamped          | SlipWall          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | Varying velocity at boundary     |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGrid-restart                  | 10 16 20     | Clamped          | Clamped           | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | Varying salt at boundary         |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | restart                          |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGridMask                      | 10 16 20     | Clamped          | Clamped           | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | Varying salt at boundary         |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | land-sea masking                 |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGridMask-1grid                | 10 16 20     | Clamped          | Clamped           | Varying salt at boundary         |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | land-sea masking                 |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGridMask-CF-Uvel              | 10 16 20     | Chapman-Flather  | Chapman-Flather   | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | Varying velocity at boundary     |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | land-sea masking                 |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGridMask-CFO-Salt             | 10 16 20     | Chapman-Flather  | Chapman-Flather   | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              | Radiation        | Radiaion          | Varying salt at boundary         |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | land-sea masking                 |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGridMask-CFO-Temp             | 10 16 20     | Chapman-Flather  | Chapman-Flather   | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              | Radiation        | Radiaion          | Varying temperature at boundary  |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | land-sea masking                 |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGridMask-CFO-Uvel             | 10 16 20     | Chapman-Flather  | Chapman-Flather   | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              | Radiation        | Radiaion          | Varying velocity at boundary     |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | land-sea masking                 |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGridMask-EWWall               | 10 16 20     | SlipWall         | Clamped           | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | Varying salt at boundary         |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | land-sea masking                 |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGridMask-EWWall-restart       | 10 16 20     | SlipWall         | Clamped           | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | Varying salt at boundary         |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | land-sea masking                 |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | restart                          |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGridMask-NSWall               | 10 16 20     | Clamped          | Slipwall          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | Varying salt at boundary         |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | land-sea masking                 |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGridMask-NSWall-restart       | 10 16 20     | Clamped          | Slipwall          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | Varying salt at boundary         |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | land-sea masking                 |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | restart                          |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGridMask-Temp                 | 10 16 20     | Clamped          | Clamped           | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | Varying temperature at boundary  |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | land-sea masking                 |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGridMask-Uvel                 | 10 16 20     | Clamped          | Clamped           | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | Varying velocity at boundary     |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | land-sea masking                 |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGridMask-Uvel-1grid           | 10 16 20     | Clamped          | Clamped           | Varying velocity at boundary     |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | land-sea masking                 |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGridMask-Uvel-EWWall          | 10 16 20     | SlipWall         | Clamped           | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | Varying salt at boundary         |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | land-sea masking                 |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGridMask-Uvel-NSWall          | 10 16 20     | Clamped          | Periodic          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | Varying salt at boundary         |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | land-sea masking                 |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGridMask-restart              | 10 16 20     | Clamped          | Clamped           | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | Varying salt at boundary         |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | land-sea masking                 |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | restart                          |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+

Based on :ref:`Particles Over Seamount<particlesseamount>`, which always includes MPI, Coriolis, and tracer particles:

+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Test                                   | nx ny nz     | xbc              | ybc               | Other                            |
+========================================+==============+==================+===================+==================================+
| ParticlesOverSeamount                  | 41 80 16     | Periodic         | Periodic          |                                  |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| ParticlesOverSeamount-restart          | 41 80 16     | Periodic         | Periodic          | restart                          |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+

Based on :ref:`Seamount<seamount>`, which always includes Coriolis and non-flat bathymetry:

+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Test                                   | nx ny nz     | xbc              | ybc               | Other                            |
+========================================+==============+==================+===================+==================================+
| Seamount-1grid-xy                      | 49 48 13     | Periodic         | Periodic          |                                  |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Seamount-xy                            | 49 48 13     | Periodic         | Periodic          | MPI                              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Seamount64-xy                          | 320 320 64   | Periodic         | Periodic          | MPI, large problem               |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+

Based on :ref:`Upwelling<upwelling>`, which always includes Coriolis and non-flat bathymetry:

+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Test                                   | nx ny nz     | xbc              | ybc               | Other                            |
+========================================+==============+==================+===================+==================================+
| Upwelling                              | 41 80 16     | Periodic         | SlipWall          | MPI                              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling-1grid                        | 41 80 16     | Periodic         | SlipWall          |                                  |
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

