# Refined where gridref_grd.nc's refinement variable asks, and its mask_rho is the land. On
# level-0 cells, land is i 16..25, j 8..14; the file asks for level 1 on i 10..21, j 3..11 and
# for level 2 on i 15..17, j 7..9, four of whose nine cells are land.
remora.max_step = 0

remora.mask_type      = "netcdf"
remora.nc_grid_file_0 = "gridref_grd.nc"

amr.max_level       = 2
amr.ref_ratio_vect  = 3 3 1 2 2 1
amr.n_error_buf     = 0
amr.blocking_factor = 1
amr.grid_eff        = 1
