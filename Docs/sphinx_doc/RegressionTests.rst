
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

The following problems are currently tested in the CI:

+----------------------+----------+----------+----------+----------+-----------------------+
| Test                 | nx ny nz | xbc      | ybc      | zbc      | Other                 |
+======================+==========+==========+==========+==========+=======================+
| Advection            | 81 81 16 | Periodic | Periodic | SlipWall |                       |
+----------------------+----------+----------+----------+----------+-----------------------+
| Advection_ML         | 80 80 16 | Periodic | Periodic | SlipWall | multilevel            |
+----------------------+----------+----------+----------+----------+-----------------------+
| DoublyPeriodic       | 41 80 16 | Periodic | Periodic | SlipWall | Coriolis              |
+----------------------+----------+----------+----------+----------+-----------------------+
| DoublyPeriodic_bathy | 41 80 16 | Periodic | Periodic | SlipWall | Coriolis              |
|                      |          |          |          |          | non-flat bathymetry   |
+----------------------+----------+----------+----------+----------+-----------------------+
| Seamount             | 49 48 13 | Periodic | Periodic | SlipWall | Coriolis              |
+----------------------+----------+----------+----------+----------+-----------------------+
| Upwelling            | 41 80 16 | Symmetry | SlipWall | SlipWall | Coriolis              |
|                      |          |          |          |          | non-flat bathymetry   |
+----------------------+----------+----------+----------+----------+-----------------------+
| Upwelling_GLS        | 41 80 16 | Periodic | SlipWall | SlipWall | Coriolis              |
|                      |          |          |          |          | non-flat bathymetry   |
|                      |          |          |          |          | GLS mixing scheme     |
+----------------------+----------+----------+----------+----------+-----------------------+

More details about the CI tests will be given below.
