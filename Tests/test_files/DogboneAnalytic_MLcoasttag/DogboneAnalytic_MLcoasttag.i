# ------------------  INPUTS TO MAIN PROGRAM  -------------------
# Asserts that a gradient criterion keyed on the mask does tag the coastline.
#
# field_name = mask is the one field exempt from the rule that a criterion sees only water:
# such a criterion is asking where the coast is, and guarding it would leave it with no
# water-water face across which the mask varies, so it would never tag anything. This is the
# documented way to refine a coastline deliberately, and it must keep working.
#
# It is also the control for DogboneAnalytic_MLcoastskip, which is this same case keyed on
# temp and asserts that nothing is tagged. Run as a pair they separate the two ways the guard
# can be wrong: without this one, a REMORAErrorTag that had quietly stopped tagging anything
# at all would still satisfy MLcoastskip.
#
# Level 1 is therefore required to exist, tracking the two coasts of the dogbone. Refining a
# coast means refining land, which the old derefine criteria would have cleared.
#
# The stationary check rides along as in DogboneAnalytic_MLmask: the exact solution is rest,
# so plt00010 must equal plt00000 with no gold file to bless. It is worth having here because
# refining across a coast is exactly where the wet-only two-way average has to hold up.
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
amr.ref_ratio_vect = 3 3  1
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

# Coast offset by one fine cell (dy_1 = 50/3) off the coarse cell face
remora.prob.mask_y_lo = 266.66666666666667
remora.prob.mask_y_hi = 483.33333333333333

# The mask is 0 or 1 exactly, so 0.5 tags precisely the cells that have a neighbor across
# the coast -- the same threshold the old internal derefine criteria used.
remora.refinement_indicators = coastline
remora.coastline.max_level = 1
remora.coastline.adjacent_difference_greater = 0.5
remora.coastline.field_name = mask

remora.coupling_type = "TwoWay"

# Mask specified on the refined level and coarsened down to level 0
remora.hires_grid_level = 1
