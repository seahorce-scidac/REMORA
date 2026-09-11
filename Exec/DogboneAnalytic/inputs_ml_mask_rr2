# ------------------  INPUTS TO MAIN PROGRAM  -------------------
# DogboneAnalytic_MLmask at an even refinement ratio, where the partially-masked blocks
# are half water rather than two thirds.
#
# Two things have to line up for one to exist. hires_grid_level resolves the coastline on the
# refined level instead of injecting it from level 0; and a static box puts refined grids over
# the coast. The box is static because this case is at rest, so a field-based indicator like
# DogboneAnalytic_MLvel's x_velocity > 0.05 would tag nothing anywhere. Offsetting
# mask_y_lo/mask_y_hi by one fine cell off a coarse face then leaves coarse rows 5 and 9
# covered by blocks that are 2 water cells out of 4.
#
# The exact solution is rest: flat bathymetry and free surface, uniform temperature and
# salinity, no initial velocity, no Coriolis. So plt00010 must equal plt00000, with no gold
# file to bless. Under the ROMS wet-only mean the partially-masked coarse cells keep
# exactly T = 10 and S = 35. Under a plain arithmetic mean the land cells, which
# advance_3d_ml zeroes every step, drag them to half of that and the run stops being
# stationary, so a lost mask weighting fails loudly.
remora.prob_name = DogboneAnalytic

remora.max_step = 10

amrex.fpe_trap_invalid=1

# PROBLEM SIZE & GEOMETRY
remora.prob_lo     =   0.       0.  -10.
remora.prob_hi     =   8400.  750.       0.

remora.n_cell           =  42 15 16

amr.blocking_factor_z = 16

amr.max_grid_size_z = 1024

remora.is_periodic = 0 0 0

remora.bc.xlo.type = "slipwall"
remora.bc.xhi.type = "slipwall"
remora.bc.ylo.type = "slipwall"
remora.bc.yhi.type = "slipwall"

# TIME STEP CONTROL
remora.fixed_dt            = 2.0 # Timestep size (seconds)
remora.ndtfast = 20

# REFINEMENT / REGRIDDING
amr.max_level       = 1       # maximum level number allowed
amr.ref_ratio_vect = 2 2  1
amr.regrid_int      = -1      # static: the grids must not move between the two plotfiles

# DIAGNOSTICS & VERBOSITY
remora.sum_interval  = 1
remora.v             = 0

# CHECKPOINT FILES
remora.check_file      = chk

# PLOTFILES
remora.plot_file     = plt
remora.plot_int      = 10
remora.plot_vars_3d  = salt temp x_velocity y_velocity z_velocity
remora.plot_vars_2d  = mask_rho   # a rho2d sidecar, for eyeballing a failure; fcompare ignores it
remora.plotfile_type = amrex

# SOLVER CHOICE
remora.tracer_horizontal_advection_scheme = "upstream3"

remora.Akt_bak = 1e-6
remora.Akv_bak = 1e-5

remora.use_coriolis  = false

remora.theta_s = 0.0
remora.theta_b = 0.0
remora.tcline = 1e16

remora.bottom_stress_type = "quadratic"
remora.rdrag2 = 3.0e-3

remora.mask_type = "analytic"

remora.init_l1ad_h =false
remora.init_l1ad_T =false

remora.init_l0int_h =true
remora.init_l0int_T =true

remora.init_ana_h = false
remora.init_ana_T = false

# PROBLEM PARAMETERS
remora.R0    = 1027.0
remora.S0    = 35.0
remora.T0    = 10.0
remora.Tcoef = 1.7e-4
remora.Scoef = 7.6e-4
remora.rho0  = 1025.0

remora.ic_type       = "analytic"

# Rest as the exact solution
remora.prob.zeta_bump = false
remora.prob.temp_west = 10.0

# Coast offset by one fine cell (dy_1 = 25) off the coarse cell face
remora.prob.mask_y_lo = 275.0
remora.prob.mask_y_hi = 475.0

# Static refinement over the coast
remora.refinement_indicators = coastbox
remora.coastbox.max_level = 1
remora.coastbox.in_box_lo = 2600. 0.
remora.coastbox.in_box_hi = 5800. 750.

remora.coupling_type = "TwoWay"

# Mask specified on the refined level and coarsened down to level 0
remora.hires_grid_level = 1
