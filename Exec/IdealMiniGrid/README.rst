This small idealized grid is used to test netCDF-provided initial and boundary
conditions. The ocean is initialized with zero velocity and a constant
temperature and salinity. Time-varying boundary conditions are then applied
for velocity, temperature, or salinity (provided by netCDF file). The default
is to use a clamped boundary condition for all quantities, but options for
Chapman-Flather and radiation conditions are available. This test also
verifies correct behavior with land-sea masking when using the ``_masked`` grid
file.

There is also a version of the problem with a file-specified surface wind
velocity (inputs file ``inputs_frc``).

To test climatology nudging, use ``inputs_clim_nudg``.

The netCDF files needed to run these tests can be found in the
`remora-data <https://github.com/seahorce-scidac/remora-data>`_
repository under the ``IdealMiniGrid`` directory.

