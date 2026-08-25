This small doubly periodic box exercises the Fennel biology model against an ocean that is
initialized and forced entirely from netCDF. The grid is 4x4 in the horizontal with 30 vertical
levels, which keeps a run short enough to compare biology tendencies term by term while still
resolving a stratified water column.

Physical state and bathymetry come from the netCDF initial and grid files. Winds, air temperature,
humidity, pressure, shortwave radiation, and rain are read from the forcing file and drive bulk
fluxes; longwave radiation is computed rather than read. Biology runs with carbon, denitrification,
and bottom-sediment fluxes enabled. Vertical mixing uses the GLS scheme, the equation of state is
nonlinear, and the bottom stress is quadratic.

The netCDF files needed to run this test can be found in the
`remora-data <https://github.com/seahorce-scidac/remora-data>`_
repository under the ``BioToy`` directory.
