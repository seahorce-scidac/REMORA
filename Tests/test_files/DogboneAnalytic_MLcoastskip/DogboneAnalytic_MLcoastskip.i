# ------------------  INPUTS TO MAIN PROGRAM  -------------------
# Asserts that a gradient criterion on a physical field does not tag the coastline.
#
# This is the case the mask guard exists for. The dogbone sits at rest with temp = 10
# throughout the water, and advance_3d_ml zeroes the tracers on land every step, so the only
# place in the whole domain where temp changes between neighbors is the coast, where it jumps
# the full 10. A criterion asking for an adjacent difference over 0.5 would tag every coastal
# cell on the strength of that jump alone, which says nothing about the flow.
#
# REMORAErrorTag takes a difference only between two water cells, so it finds no difference
# anywhere and tags nothing: the run must stay single-level. That is what check_max_level.sh
# asserts, and it is the whole assertion -- fcompare cannot see it, because a run that wrongly
# refined the coast would still produce a perfectly self-consistent plotfile.
#
# Its companion DogboneAnalytic_MLcoasttag is the same case keyed on the mask instead, which
# is exempt from the guard and does build level 1 along the coast. The pair is what
# distinguishes "correctly declined to tag" from "quietly stopped tagging anything".
#
# The stationary check rides along for free, as in DogboneAnalytic_MLmask: the exact solution
# is rest, so plt00010 must equal plt00000 with no gold file to bless.
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

# Threshold well under the 10 the coast jumps by, so nothing but the guard can keep this
# criterion from tagging there.
remora.refinement_indicators = coastgrad
remora.coastgrad.max_level = 1
remora.coastgrad.adjacent_difference_greater = 0.5
remora.coastgrad.field_name = temp

remora.coupling_type = "TwoWay"

# Mask specified on the refined level and coarsened down to level 0
remora.hires_grid_level = 1
