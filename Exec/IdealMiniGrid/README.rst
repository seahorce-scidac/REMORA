This small idealized grid is used to test netCDF-provided initial and boundary
conditions. The ocean is initialized with zero velocity and a constant
temperature and salinity. Time-varying boundary conditions are then applied
for velocity, temperature, or salinity (provided by netCDF file). The default
is to use a clamped boundary condition for all quantities, but options for
Chapman-Flather and radiation conditions are available. This test also
verifies correct behavior with land-sea masking when using the ``_masked`` grid
file.

There is also a version of the problem with a file-specified surface wind
velocity (inputs file ``inputs_frc``). ``inputs_frc_latlon`` gives the winds and
air temperature on a coarse lat-lon grid instead, over a rotated copy of the
grid, so that they are interpolated onto the model grid and the winds rotated
onto its axes; its header gives the commands that make its two files with
``Tests/tools/make_latlon_forcing_test_data.py`` and check the result.

To test climatology nudging, use ``inputs_clim_nudg``.

``inputs_synth_hires`` reads the grid and initial state on level 1 from synthetic files written
by ``Tests/tools/make_hires_test_data.py``. ``inputs_synth_hires_refine`` places the refined
grids from the ``refinement`` variable of the high-resolution grid file, on two levels; run
the generator with ``--with-refinement`` first, and compare the grids the run prints with the
ones in the file's header.

The netCDF files needed to run these tests can be found in the
`remora-data <https://github.com/seahorce-scidac/remora-data>`_
repository under the ``IdealMiniGrid`` directory.

