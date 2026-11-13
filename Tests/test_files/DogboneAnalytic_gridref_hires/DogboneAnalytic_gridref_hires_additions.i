# As DogboneAnalytic_gridref, with the mask and the refinement variable given on level 1
# (gridref_hires_grd.nc, ratio 3). Land is the fine cells under level-0 cells i 16..25,
# j 8..14. The file asks for level 1 on fine cells i 31..64, j 10..34, which start and end
# partway into level-0 cells 10..21 and 3..11.
remora.max_step = 0

remora.mask_type          = "netcdf"
remora.hires_grid_level   = 1
remora.nc_grid_file_hires = "gridref_hires_grd.nc"

amr.max_level       = 1
amr.ref_ratio_vect  = 3 3 1
amr.n_error_buf     = 0
amr.blocking_factor = 1
amr.grid_eff        = 1
