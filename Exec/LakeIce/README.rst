LakeIce
=======

The myROMS LAKE_ICE test (myroms/roms_test ``lake_ice``) with REMORA's port of
the myROMS sea-ice model: a 200 km x 100 km lake, 20 to 200 m deep, under the
analytic seasonal atmosphere of the ROMS case (``Prob/REMORA_InitAnalyticSurfaceVar_LakeIce.H``).
These are the inputs of the REMORA-ROMS comparison described in the Sea Ice
chapter of the documentation.

The inputs start from the day-60 state of a ROMS spin-up, when a quarter of the
lake is ice covered, and need two NetCDF files that are not in the repository:

* ``roms_lake_ice_grd.nc``: the grid of the ROMS case, from roms_test
  ``lake_ice/Data`` (a git-lfs file).
* ``lake_ice_day60_ini_remora.nc``: the initial state. Run ROMS LAKE_ICE (built
  with ``ICE_MODEL BULK_FLUXES GLS_MIXING`` in place of ``LMD_MIXING``) for 60
  days from its ice-free initial file, and take the last record of the restart
  file. Drop ``Uevp`` and ``Vevp`` (ROMS cannot read its own restart with them)
  and replace the fill values on land with zero.

REMORA must be built with NetCDF (``REMORA_ENABLE_PNETCDF``). The inputs set
``remora.bulk_flux_scheme = coare35`` and ``remora.g = 9.81`` to match myROMS.
