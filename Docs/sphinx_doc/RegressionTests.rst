
 .. _RegressionTests:

Regression Tests
================

There are currently 8 tests which are run as part of every PR.
The CI tests use cmake and are based on the version
of AMReX in the REMORA submodule.

In addition there are nightly tests that use GNUMake and use the current
development branch of AMReX.

Results from the nightly CPU tests can be found here: `CPU tests`_

Results from the nightly GPU tests can be found here: `GPU tests`_

.. _`CPU tests`: https://ccse.lbl.gov/pub/RegressionTesting1/REMORA

.. _`GPU tests`: https://ccse.lbl.gov/pub/GpuRegressionTesting/REMORA

The following problems are currently tested in the CI. More details about the CI tests are given below.

+----------------------+----------+----------+----------+-----------------------+
| Test                 | nx ny nz | xbc      | ybc      | Other                 |
+======================+==========+==========+==========+=======================+
| Advection            | 81 81 16 | Periodic | Periodic |                       |
+----------------------+----------+----------+----------+-----------------------+
| Advection_ML         | 80 80 16 | Periodic | Periodic | multilevel            |
+----------------------+----------+----------+----------+-----------------------+
| Channel_Test         | 20 60 50 | Periodic | SlipWall | Coriolis              |
|                      |          |          |          | GLS mixing scheme     |
|                      |          |          |          | non-flat bathymetry   |
+----------------------+----------+----------+----------+-----------------------+
| DoubleGyre           | 54 108 4 | SlipWall | SlipWall | Coriolis              |
+----------------------+----------+----------+----------+-----------------------+
| DoublyPeriodic       | 41 80 16 | Periodic | Periodic | Coriolis              |
+----------------------+----------+----------+----------+-----------------------+
| DoublyPeriodic_bathy | 41 80 16 | Periodic | Periodic | Coriolis              |
|                      |          |          |          | non-flat bathymetry   |
+----------------------+----------+----------+----------+-----------------------+
| Seamount             | 49 48 13 | Periodic | Periodic | Coriolis              |
+----------------------+----------+----------+----------+-----------------------+
| Upwelling            | 41 80 16 | Periodic | SlipWall | Coriolis              |
|                      |          |          |          | non-flat bathymetry   |
+----------------------+----------+----------+----------+-----------------------+
| Upwelling_GLS        | 41 80 16 | Periodic | SlipWall | Coriolis              |
|                      |          |          |          | non-flat bathymetry   |
|                      |          |          |          | GLS mixing scheme     |
+----------------------+----------+----------+----------+-----------------------+

And the following are currently tested nighly on CPU

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
|                                        |              |                  |                   | multilevel                       |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Advection-xy                           | 81 81 16     | Periodic         | Periodic          | MPI                              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Advection-xy-restart                   | 81 81 16     | Periodic         | Periodic          | MPI                              |
|                                        |              |                  |                   | restart                          |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Advection-xy-restart-ML                | 80 80 16     | Periodic         | Periodic          | MPI                              |
|                                        |              |                  |                   | restart                          |
|                                        |              |                  |                   | multilevel                       |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Advection-U3-xy                        | 81 81 16     | Periodic         | Periodic          | MPI                              |
|                                        |              |                  |                   | advection: upstream 3rd order    |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| ChannelTest                            | 20 60 50     | Periodic         | SlipWall          | MPI                              |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | GLS mixing scheme                |
|                                        |              |                  |                   | non-flat bathymetry              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| ChannelTest-OMP                        | 20 60 50     | Periodic         | SlipWall          | MPI + OpenMP                     |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | GLS mixing scheme                |
|                                        |              |                  |                   | non-flat bathymetry              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| ChannelTest-OMP-1grid-xy               | 20 60 50     | Periodic         | SlipWall          | OpenMP                           |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | GLS mixing scheme                |
|                                        |              |                  |                   | non-flat bathymetry              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| ChannelTest-xy-restart                 | 20 60 50     | Periodic         | SlipWall          | MPI + OpenMP                     |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | GLS mixing scheme                |
|                                        |              |                  |                   | non-flat bathymetry              |
|                                        |              |                  |                   | restart                          |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| ChannelTestOrlanski                    | 20 60 50     | Radiation        | Radiation         | MPI                              |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | GLS mixing scheme                |
|                                        |              |                  |                   | non-flat bathymetry              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| ChannelTestOrlanski-OMP                | 20 60 50     | Radiation        | Radiation         | MPI + OpenMP                     |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | GLS mixing scheme                |
|                                        |              |                  |                   | non-flat bathymetry              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| ChannelTestOrlanski-OMP-1grid-xy       | 20 60 50     | Radiation        | Radiation         | OpenMP                           |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | GLS mixing scheme                |
|                                        |              |                  |                   | non-flat bathymetry              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| ChannelTestOrlanski-xy-restart         | 20 60 50     | Radiation        | Radiation         | MPI + OpenMP                     |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | GLS mixing scheme                |
|                                        |              |                  |                   | non-flat bathymetry              |
|                                        |              |                  |                   | restart                          |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DoubleGyre                             | 54 108 4     | SlipWall         | SlipWall          | MPI                              |
|                                        |              |                  |                   | Coriolis                         |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DoubleGyre-OMP                         | 54 108 4     | SlipWall         | SlipWall          | MPI + OpenMP                     |
|                                        |              |                  |                   | Coriolis                         |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DoubleGyre-OMP-1grid-xy                | 54 108 4     | SlipWall         | SlipWall          | OpenMP                           |
|                                        |              |                  |                   | Coriolis                         |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DoubleGyre-OMP-xy-restart              | 54 108 4     | SlipWall         | SlipWall          | MPI + OpenMP                     |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | restart                          |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DoublyPeriodic-1grid-xy                | 41 80 16     | Periodic         | Periodic          | Coriolis                         |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DoublyPeriodic-1grid-xy-bathy          | 41 80 16     | Periodic         | Periodic          | Coriolis                         |
|                                        |              |                  |                   | non-flat bathymetry              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DoublyPeriodic-NETCDF-build            | 41 80 16     | N/A              | N/A               | Coriolis                         |
|                                        |              |                  |                   | Build w/PnetCDF                  |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DoublyPeriodic-OMP-1grid-xy            | 41 80 16     | Periodic         | Periodic          | OpenMP                           |
|                                        |              |                  |                   | Coriolis                         |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DoublyPeriodic-OMP-xy                  | 41 80 16     | Periodic         | Periodic          | MPI + OpenMP                     |
|                                        |              |                  |                   | Coriolis                         |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DoublyPeriodic-OMP-xy-bathy            | 41 80 16     | Periodic         | Periodic          | MPI + OpenMP                     |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | non-flat bathymetry              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DoublyPeriodic-xy                      | 41 80 16     | Periodic         | Periodic          | MPI                              |
|                                        |              |                  |                   | Coriolis                         |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DoublyPeriodic-xy-bathy                | 41 80 16     | Periodic         | Periodic          | MPI                              |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | non-flat bathymetry              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DoublyPeriodic-xy-restart              | 41 80 16     | Periodic         | Periodic          | MPI                              |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | restart                          |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DoublyPeriodic64-OMP-xy                | 328 320 64   | Periodic         | Periodic          | MPI + OpenMP, large problem      |
|                                        |              |                  |                   | Coriolis                         |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DoublyPeriodic64-OMP-xy-bathy          | 328 320 64   | Periodic         | Periodic          | MPI + OpenMP, large problem      |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | non-flat bathymetry              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DoublyPeriodic64-xy                    | 328 320 64   | Periodic         | Periodic          | MPI, large problem               |
|                                        |              |                  |                   | Coriolis                         |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DoublyPeriodicC4-xy                    | 41 80 16     | Periodic         | Periodic          | MPI                              |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | advection: centered 4th order    |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGrid                          | 10 16 20     | Clamped          | Clamped           | MPI                              |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | PnetCDF                          |
|                                        |              |                  |                   | Varying salt at boundary         |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGrid-1grid                    | 10 16 20     | Clamped          | Clamped           | Coriolis                         |
|                                        |              |                  |                   | PnetCDF                          |
|                                        |              |                  |                   | Varying salt at boundary         |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGrid-CF-Uvel-OMP              | 10 16 20     | Chapman-Flather  | Chapman-Flather   | MPI + OpenMP                     |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | PnetCDF                          |
|                                        |              |                  |                   | Varying velocity at boundary     |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGrid-CFO-Salt-OMP             | 10 16 20     | Chapman-Flather  | Chapman-Flather   | MPI + OpenMP                     |
|                                        |              | Radiation        | Radiation         | Coriolis                         |
|                                        |              |                  |                   | PnetCDF                          |
|                                        |              |                  |                   | Varying salt at boundary         |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGrid-CFO-Temp-OMP             | 10 16 20     | Chapman-Flather  | Chapman-Flather   | MPI + OpenMP                     |
|                                        |              | Radiation        | Radiation         | Coriolis                         |
|                                        |              |                  |                   | PnetCDF                          |
|                                        |              |                  |                   | Varying temperature at boundary  |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGrid-CFO-Uvel-OMP             | 10 16 20     | Chapman-Flather  | Chapman-Flather   | MPI + OpenMP                     |
|                                        |              | Radiation        | Radiation         | Coriolis                         |
|                                        |              |                  |                   | PnetCDF                          |
|                                        |              |                  |                   | Varying velocity at boundary     |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGrid-EWWall-OMP               | 10 16 20     | SlipWall         | Clamped           | MPI + OpenMP                     |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | PnetCDF                          |
|                                        |              |                  |                   | Varying salt at boundary         |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGrid-EWWall-restart           | 10 16 20     | SlipWall         | Clamped           | MPI                              |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | PnetCDF                          |
|                                        |              |                  |                   | Varying salt at boundary         |
|                                        |              |                  |                   | restart                          |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGrid-NSWall-OMP               | 10 16 20     | Clamped          | SlipWall          | MPI + OpenMP                     |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | PnetCDF                          |
|                                        |              |                  |                   | Varying salt at boundary         |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGrid-NSWall-restart           | 10 16 20     | Clamped          | Slipwall          | MPI                              |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | PnetCDF                          |
|                                        |              |                  |                   | Varying salt at boundary         |
|                                        |              |                  |                   | restart                          |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGrid-OMP                      | 10 16 20     | Clamped          | Clamped           | MPI + OpenMP                     |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | PnetCDF                          |
|                                        |              |                  |                   | Varying salt at boundary         |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGrid-OMP-1grid                | 10 16 20     | Clamped          | Clamped           | OpenMP                           |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | PnetCDF                          |
|                                        |              |                  |                   | Varying salt at boundary         |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGrid-Temp                     | 10 16 20     | Clamped          | Clamped           | MPI                              |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | PnetCDF                          |
|                                        |              |                  |                   | Varying temperature at boundary  |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGrid-Uvel                     | 10 16 20     | Clamped          | Clamped           | MPI                              |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | PnetCDF                          |
|                                        |              |                  |                   | Varying velocity at boundary     |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGrid-Uvel-1grid               | 10 16 20     | Clamped          | Clamped           | Coriolis                         |
|                                        |              |                  |                   | PnetCDF                          |
|                                        |              |                  |                   | Varying velocity at boundary     |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGrid-Uvel-EWWall-OMP          | 10 16 20     | SlipWall         | Clamped           | MPI + OpenMP                     |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | PnetCDF                          |
|                                        |              |                  |                   | Varying velocity at boundary     |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGrid-Uvel-NSWall-OMP          | 10 16 20     | Clamped          | SlipWall          | MPI + OpenMP                     |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | PnetCDF                          |
|                                        |              |                  |                   | Varying velocity at boundary     |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGrid-Uvel-OMP                 | 10 16 20     | Clamped          | Clamped           | MPI + OpenMP                     |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | PnetCDF                          |
|                                        |              |                  |                   | Varying velocity at boundary     |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGrid-Uvel-OMP-1grid           | 10 16 20     | Clamped          | Clamped           | OpenMP                           |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | PnetCDF                          |
|                                        |              |                  |                   | Varying velocity at boundary     |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGrid-restart                  | 10 16 20     | Clamped          | Clamped           | MPI                              |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | PnetCDF                          |
|                                        |              |                  |                   | Varying salt at boundary         |
|                                        |              |                  |                   | restart                          |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGridMask                      | 10 16 20     | Clamped          | Clamped           | MPI                              |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | PnetCDF                          |
|                                        |              |                  |                   | Varying salt at boundary         |
|                                        |              |                  |                   | land-sea masking                 |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGridMask-1grid                | 10 16 20     | Clamped          | Clamped           | Coriolis                         |
|                                        |              |                  |                   | PnetCDF                          |
|                                        |              |                  |                   | Varying salt at boundary         |
|                                        |              |                  |                   | land-sea masking                 |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGridMask-CF-Uvel-OMP          | 10 16 20     | Chapman-Flather  | Chapman-Flather   | MPI + OpenMP                     |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | PnetCDF                          |
|                                        |              |                  |                   | Varying velocity at boundary     |
|                                        |              |                  |                   | land-sea masking                 |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGridMask-CFO-Salt-OMP         | 10 16 20     | Chapman-Flather  | Chapman-Flather   | MPI + OpenMP                     |
|                                        |              | Radiation        | Radiaion          | Coriolis                         |
|                                        |              |                  |                   | PnetCDF                          |
|                                        |              |                  |                   | Varying salt at boundary         |
|                                        |              |                  |                   | land-sea masking                 |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGridMask-CFO-Temp-OMP         | 10 16 20     | Chapman-Flather  | Chapman-Flather   | MPI + OpenMP                     |
|                                        |              | Radiation        | Radiaion          | Coriolis                         |
|                                        |              |                  |                   | PnetCDF                          |
|                                        |              |                  |                   | Varying temperature at boundary  |
|                                        |              |                  |                   | land-sea masking                 |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGridMask-CFO-Uvel-OMP         | 10 16 20     | Chapman-Flather  | Chapman-Flather   | MPI + OpenMP                     |
|                                        |              | Radiation        | Radiaion          | Coriolis                         |
|                                        |              |                  |                   | PnetCDF                          |
|                                        |              |                  |                   | Varying velocity at boundary     |
|                                        |              |                  |                   | land-sea masking                 |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGridMask-EWWall-OMP           | 10 16 20     | SlipWall         | Clamped           | MPI + OpenMP                     |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | PnetCDF                          |
|                                        |              |                  |                   | Varying salt at boundary         |
|                                        |              |                  |                   | land-sea masking                 |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGridMask-OMP                  | 10 16 20     | Clamped          | Clamped           | MPI + OpenMP                     |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | PnetCDF                          |
|                                        |              |                  |                   | Varying salt at boundary         |
|                                        |              |                  |                   | land-sea masking                 |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGridMask-OMP-1grid            | 10 16 20     | Clamped          | Clamped           | OpenMP                           |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | PnetCDF                          |
|                                        |              |                  |                   | Varying salt at boundary         |
|                                        |              |                  |                   | land-sea masking                 |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGridMask-Temp                 | 10 16 20     | Clamped          | Clamped           | MPI                              |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | PnetCDF                          |
|                                        |              |                  |                   | Varying temperature at boundary  |
|                                        |              |                  |                   | land-sea masking                 |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGridMask-Uvel                 | 10 16 20     | Clamped          | Clamped           | MPI                              |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | PnetCDF                          |
|                                        |              |                  |                   | Varying velocity at boundary     |
|                                        |              |                  |                   | land-sea masking                 |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGridMask-Uvel-1grid           | 10 16 20     | Clamped          | Clamped           | Coriolis                         |
|                                        |              |                  |                   | PnetCDF                          |
|                                        |              |                  |                   | Varying velocity at boundary     |
|                                        |              |                  |                   | land-sea masking                 |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGridMask-Uvel-EWWall-OMP      | 10 16 20     | SlipWall         | Clamped           | MPI + OpenMP                     |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | PnetCDF                          |
|                                        |              |                  |                   | Varying velocity at boundary     |
|                                        |              |                  |                   | land-sea masking                 |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGridMask-Uvel-NSWall-OMP      | 10 16 20     | Clamped          | SlipWall          | MPI + OpenMP                     |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | PnetCDF                          |
|                                        |              |                  |                   | Varying velocity at boundary     |
|                                        |              |                  |                   | land-sea masking                 |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGridMask-Uvel-OMP             | 10 16 20     | Clamped          | Clamped           | MPI + OpenMP                     |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | PnetCDF                          |
|                                        |              |                  |                   | Varying velocity at boundary     |
|                                        |              |                  |                   | land-sea masking                 |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGridMask-Uvel-OMP-1grid       | 10 16 20     | Clamped          | Clamped           | OpenMP                           |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | PnetCDF                          |
|                                        |              |                  |                   | Varying velocity at boundary     |
|                                        |              |                  |                   | land-sea masking                 |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGridMask-restart              | 10 16 20     | Clamped          | Clamped           | MPI                              |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | PnetCDF                          |
|                                        |              |                  |                   | Varying salt at boundary         |
|                                        |              |                  |                   | land-sea masking                 |
|                                        |              |                  |                   | restart                          |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| ParticlesOverSeamount                  | 41 80 16     | Periodic         | Periodic          | MPI                              |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | Tracer particles                 |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| ParticlesOverSeamount-restart          | 41 80 16     | Periodic         | Periodic          | MPI                              |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | Tracer particles                 |
|                                        |              |                  |                   | restart                          |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Seamount-1grid-xy                      | 49 48 13     | Periodic         | Periodic          | Coriolis                         |
|                                        |              |                  |                   | non-flat bathymetry              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Seamount-OMP-xy                        | 49 48 13     | Periodic         | Periodic          | MPI + OpenMP                     |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | non-flat bathymetry              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Seamount-xy                            | 49 48 13     | Periodic         | Periodic          | MPI                              |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | non-flat bathymetry              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Seamount-xy-restart                    | 49 48 13     | Periodic         | Periodic          | MPI                              |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | non-flat bathymetry              |
|                                        |              |                  |                   | restart                          |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Seamount64-OMP-xy                      | 320 320 64   | Periodic         | Periodic          | MPI + OpenMP, large problem      |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | non-flat bathymetry              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling                              | 41 80 16     | Periodic         | SlipWall          | MPI                              |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | non-flat bathymetry              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling-1grid                        | 41 80 16     | Periodic         | SlipWall          | Coriolis                         |
|                                        |              |                  |                   | non-flat bathymetry              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling-OMP                          | 41 80 16     | Periodic         | SlipWall          | MPI + OpenMP                     |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | non-flat bathymetry              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling-OMP-1grid                    | 41 80 16     | Periodic         | SlipWall          | OpenMP                           |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | non-flat bathymetry              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling-restart                      | 41 80 16     | Periodic         | SlipWall          | MPI                              |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | non-flat bathymetry              |
|                                        |              |                  |                   | restart                          |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling-x                            | 41 80 16     | SlipWall         | Periodic          | MPI                              |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | non-flat bathymetry              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling-x-1grid                      | 41 80 16     | SlipWall         | Periodic          | Coriolis                         |
|                                        |              |                  |                   | non-flat bathymetry              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling-x-OMP                        | 41 80 16     | SlipWall         | Periodic          | MPI + OpenMP                     |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | non-flat bathymetry              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling64-OMP                        | 328 320 64   | SlipWall         | Periodic          | MPI + OpenMP, large problem      |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | non-flat bathymetry              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling64-OMP                        | 328 320 64   | SlipWall         | Periodic          | MPI + OpenMP, large problem      |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | non-flat bathymetry              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| UpwellingC4                            | 41 80 16     | Periodic         | SlipWall          | MPI                              |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | non-flat bathymetry              |
|                                        |              |                  |                   | advection: centered 4th order    |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling_GLS                          | 41 80 16     | Periodic         | SlipWall          | MPI                              |
|                                        |              |                  |                   | GLS mixing scheme                |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | non-flat bathymetry              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling_GLS-restart                  | 41 80 16     | Periodic         | SlipWall          | MPI                              |
|                                        |              |                  |                   | GLS mixing scheme                |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | non-flat bathymetry              |
|                                        |              |                  |                   | restart                          |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling_GLS_Canuto_A                 | 41 80 16     | Periodic         | SlipWall          | MPI                              |
|                                        |              |                  |                   | GLS mixing scheme                |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | non-flat bathymetry              |
|                                        |              |                  |                   | Canuto A stability               |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling_GLS_Canuto_B                 | 41 80 16     | Periodic         | SlipWall          | MPI                              |
|                                        |              |                  |                   | GLS mixing scheme                |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | non-flat bathymetry              |
|                                        |              |                  |                   | Canuto B stability               |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+

And the following are currently tested nighly on GPU. All are compiled and run with CUDA

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
|                                        |              |                  |                   | multilevel                       |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Advection-xy-restart                   | 81 81 16     | Periodic         | Periodic          | MPI                              |
|                                        |              |                  |                   | restart                          |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Advection64-xy                         | 328 320 64   | Periodic         | Periodic          | MPI, large problem               |
|                                        |              |                  |                   | restart                          |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Advection-U3-xy                        | 81 81 16     | Periodic         | Periodic          | MPI                              |
|                                        |              |                  |                   | advection: upstream 3rd order    |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| ChannelTest-1grid-xy                   | 20 60 50     | Periodic         | SlipWall          | Coriolis                         |
|                                        |              |                  |                   | GLS mixing scheme                |
|                                        |              |                  |                   | non-flat bathymetry              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| ChannelTest-xy                         | 20 60 50     | Periodic         | SlipWall          | MPI                              |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | GLS mixing scheme                |
|                                        |              |                  |                   | non-flat bathymetry              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| ChannelTest-xy-restart                 | 20 60 50     | Periodic         | SlipWall          | MPI                              |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | GLS mixing scheme                |
|                                        |              |                  |                   | non-flat bathymetry              |
|                                        |              |                  |                   | restart                          |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DoubleGyre-1grid-xy                    | 54 108 4     | SlipWall         | SlipWall          | Coriolis                         |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DoubleGyre-xy                          | 54 108 4     | SlipWall         | SlipWall          | MPI                              |
|                                        |              |                  |                   | Coriolis                         |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DoubleGyre-xy-restart                  | 54 108 4     | SlipWall         | SlipWall          | MPI                              |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | restart                          |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DoublyPeriodic-1grid-xy                | 41 80 16     | Periodic         | Periodic          | Coriolis                         |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DoublyPeriodic-xy                      | 41 80 16     | Periodic         | Periodic          | MPI                              |
|                                        |              |                  |                   | Coriolis                         |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DoublyPeriodic-xy-bathy                | 41 80 16     | Periodic         | Periodic          | MPI                              |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | non-flat bathymetry              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DoublyPeriodic-xy-restart              | 41 80 16     | Periodic         | Periodic          | MPI                              |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | restart                          |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DoublyPeriodic64-xy                    | 328 320 64   | Periodic         | Periodic          | MPI, large problem               |
|                                        |              |                  |                   | Coriolis                         |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DoublyPeriodic64-xy-bathy              | 328 320 64   | Periodic         | Periodic          | MPI, large problem               |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | non-flat bathymetry              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DoublyPeriodicC4-xy                    | 41 80 16     | Periodic         | Periodic          | MPI                              |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | advection: centered 4th order    |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGrid                          | 10 16 20     | Clamped          | Clamped           | MPI                              |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | PnetCDF                          |
|                                        |              |                  |                   | Varying salt at boundary         |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGrid-1grid                    | 10 16 20     | Clamped          | Clamped           | Coriolis                         |
|                                        |              |                  |                   | PnetCDF                          |
|                                        |              |                  |                   | Varying salt at boundary         |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGrid-CF-Uvel                  | 10 16 20     | Chapman-Flather  | Chapman-Flather   | MPI                              |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | PnetCDF                          |
|                                        |              |                  |                   | Varying velocity at boundary     |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGrid-CFO-Salt                 | 10 16 20     | Chapman-Flather  | Chapman-Flather   | MPI                              |
|                                        |              | Radiation        | Radiation         | Coriolis                         |
|                                        |              |                  |                   | PnetCDF                          |
|                                        |              |                  |                   | Varying salt at boundary         |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGrid-CFO-Temp                 | 10 16 20     | Chapman-Flather  | Chapman-Flather   | MPI                              |
|                                        |              | Radiation        | Radiation         | Coriolis                         |
|                                        |              |                  |                   | PnetCDF                          |
|                                        |              |                  |                   | Varying temperature at boundary  |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGrid-CFO-Uvel                 | 10 16 20     | Chapman-Flather  | Chapman-Flather   | MPI                              |
|                                        |              | Radiation        | Radiation         | Coriolis                         |
|                                        |              |                  |                   | PnetCDF                          |
|                                        |              |                  |                   | Varying velocity at boundary     |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGrid-EWWall                   | 10 16 20     | SlipWall         | Clamped           | MPI                              |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | PnetCDF                          |
|                                        |              |                  |                   | Varying salt at boundary         |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGrid-NSWall                   | 10 16 20     | Clamped          | SlipWall          | MPI                              |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | PnetCDF                          |
|                                        |              |                  |                   | Varying salt at boundary         |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGrid-Temp                     | 10 16 20     | Clamped          | Clamped           | MPI                              |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | PnetCDF                          |
|                                        |              |                  |                   | Varying temperature at boundary  |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGrid-Uvel                     | 10 16 20     | Clamped          | Clamped           | MPI                              |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | PnetCDF                          |
|                                        |              |                  |                   | Varying velocity at boundary     |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGrid-Uvel-1grid               | 10 16 20     | Clamped          | Clamped           | Coriolis                         |
|                                        |              |                  |                   | PnetCDF                          |
|                                        |              |                  |                   | Varying velocity at boundary     |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGrid-Uvel-EWWall              | 10 16 20     | SlipWall         | Clamped           | MPI                              |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | PnetCDF                          |
|                                        |              |                  |                   | Varying velocity at boundary     |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGrid-Uvel-NSWall              | 10 16 20     | Clamped          | SlipWall          | MPI                              |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | PnetCDF                          |
|                                        |              |                  |                   | Varying velocity at boundary     |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGrid-restart                  | 10 16 20     | Clamped          | Clamped           | MPI                              |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | PnetCDF                          |
|                                        |              |                  |                   | Varying salt at boundary         |
|                                        |              |                  |                   | restart                          |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGridMask                      | 10 16 20     | Clamped          | Clamped           | MPI                              |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | PnetCDF                          |
|                                        |              |                  |                   | Varying salt at boundary         |
|                                        |              |                  |                   | land-sea masking                 |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGridMask-1grid                | 10 16 20     | Clamped          | Clamped           | Coriolis                         |
|                                        |              |                  |                   | PnetCDF                          |
|                                        |              |                  |                   | Varying salt at boundary         |
|                                        |              |                  |                   | land-sea masking                 |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGridMask-CF-Uvel              | 10 16 20     | Chapman-Flather  | Chapman-Flather   | MPI                              |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | PnetCDF                          |
|                                        |              |                  |                   | Varying velocity at boundary     |
|                                        |              |                  |                   | land-sea masking                 |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGridMask-CFO-Salt             | 10 16 20     | Chapman-Flather  | Chapman-Flather   | MPI                              |
|                                        |              | Radiation        | Radiaion          | Coriolis                         |
|                                        |              |                  |                   | PnetCDF                          |
|                                        |              |                  |                   | Varying salt at boundary         |
|                                        |              |                  |                   | land-sea masking                 |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGridMask-CFO-Temp             | 10 16 20     | Chapman-Flather  | Chapman-Flather   | MPI                              |
|                                        |              | Radiation        | Radiaion          | Coriolis                         |
|                                        |              |                  |                   | PnetCDF                          |
|                                        |              |                  |                   | Varying temperature at boundary  |
|                                        |              |                  |                   | land-sea masking                 |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGridMask-CFO-Uvel             | 10 16 20     | Chapman-Flather  | Chapman-Flather   | MPI                              |
|                                        |              | Radiation        | Radiaion          | Coriolis                         |
|                                        |              |                  |                   | PnetCDF                          |
|                                        |              |                  |                   | Varying velocity at boundary     |
|                                        |              |                  |                   | land-sea masking                 |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGridMask-EWWall               | 10 16 20     | SlipWall         | Clamped           | MPI                              |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | PnetCDF                          |
|                                        |              |                  |                   | Varying salt at boundary         |
|                                        |              |                  |                   | land-sea masking                 |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGridMask-EWWall-restart       | 10 16 20     | SlipWall         | Clamped           | MPI                              |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | PnetCDF                          |
|                                        |              |                  |                   | Varying salt at boundary         |
|                                        |              |                  |                   | land-sea masking                 |
|                                        |              |                  |                   | restart                          |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGridMask-NSWall               | 10 16 20     | Clamped          | Slipwall          | MPI                              |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | PnetCDF                          |
|                                        |              |                  |                   | Varying salt at boundary         |
|                                        |              |                  |                   | land-sea masking                 |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGridMask-NSWall-restart       | 10 16 20     | Clamped          | Slipwall          | MPI                              |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | PnetCDF                          |
|                                        |              |                  |                   | Varying salt at boundary         |
|                                        |              |                  |                   | land-sea masking                 |
|                                        |              |                  |                   | restart                          |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGridMask-Temp                 | 10 16 20     | Clamped          | Clamped           | MPI                              |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | PnetCDF                          |
|                                        |              |                  |                   | Varying temperature at boundary  |
|                                        |              |                  |                   | land-sea masking                 |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGridMask-Uvel                 | 10 16 20     | Clamped          | Clamped           | MPI                              |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | PnetCDF                          |
|                                        |              |                  |                   | Varying velocity at boundary     |
|                                        |              |                  |                   | land-sea masking                 |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGridMask-Uvel-1grid           | 10 16 20     | Clamped          | Clamped           | Coriolis                         |
|                                        |              |                  |                   | PnetCDF                          |
|                                        |              |                  |                   | Varying velocity at boundary     |
|                                        |              |                  |                   | land-sea masking                 |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGridMask-Uvel-EWWall          | 10 16 20     | SlipWall         | Clamped           | MPI                              |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | PnetCDF                          |
|                                        |              |                  |                   | Varying salt at boundary         |
|                                        |              |                  |                   | land-sea masking                 |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGridMask-Uvel-NSWall          | 10 16 20     | Clamped          | Periodic          | MPI                              |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | PnetCDF                          |
|                                        |              |                  |                   | Varying salt at boundary         |
|                                        |              |                  |                   | land-sea masking                 |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| IdealMiniGridMask-restart              | 10 16 20     | Clamped          | Clamped           | MPI                              |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | PnetCDF                          |
|                                        |              |                  |                   | Varying salt at boundary         |
|                                        |              |                  |                   | land-sea masking                 |
|                                        |              |                  |                   | restart                          |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| ParticlesOverSeamount                  | 41 80 16     | Periodic         | Periodic          | MPI                              |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | Tracer particles                 |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| ParticlesOverSeamount-restart          | 41 80 16     | Periodic         | Periodic          | MPI                              |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | Tracer particles                 |
|                                        |              |                  |                   | restart                          |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Seamount-1grid-xy                      | 49 48 13     | Periodic         | Periodic          | Coriolis                         |
|                                        |              |                  |                   | non-flat bathymetry              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Seamount-xy                            | 49 48 13     | Periodic         | Periodic          | MPI                              |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | non-flat bathymetry              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Seamount64-xy                          | 320 320 64   | Periodic         | Periodic          | MPI, large problem               |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | non-flat bathymetry              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling                              | 41 80 16     | Periodic         | SlipWall          | MPI                              |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | non-flat bathymetry              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling-1grid                        | 41 80 16     | Periodic         | SlipWall          | Coriolis                         |
|                                        |              |                  |                   | non-flat bathymetry              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling-x                            | 41 80 16     | SlipWall         | Periodic          | MPI                              |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | non-flat bathymetry              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling-x-1grid                      | 41 80 16     | SlipWall         | Periodic          | Coriolis                         |
|                                        |              |                  |                   | non-flat bathymetry              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling64                            | 328 320 64   | SlipWall         | Periodic          | MPI, large problem               |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | non-flat bathymetry              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| UpwellingC4                            | 41 80 16     | Periodic         | SlipWall          | MPI                              |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | non-flat bathymetry              |
|                                        |              |                  |                   | advection: centered 4th order    |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling_GLS                          | 41 80 16     | Periodic         | SlipWall          | MPI                              |
|                                        |              |                  |                   | GLS mixing scheme                |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | non-flat bathymetry              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling_GLS-restart                  | 41 80 16     | Periodic         | SlipWall          | MPI                              |
|                                        |              |                  |                   | GLS mixing scheme                |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | non-flat bathymetry              |
|                                        |              |                  |                   | restart                          |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling_GLS_Canuto_A                 | 41 80 16     | Periodic         | SlipWall          | MPI                              |
|                                        |              |                  |                   | GLS mixing scheme                |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | non-flat bathymetry              |
|                                        |              |                  |                   | Canuto A stability               |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling_GLS_Canuto_B                 | 41 80 16     | Periodic         | SlipWall          | MPI                              |
|                                        |              |                  |                   | GLS mixing scheme                |
|                                        |              |                  |                   | Coriolis                         |
|                                        |              |                  |                   | non-flat bathymetry              |
|                                        |              |                  |                   | Canuto B stability               |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+

