# ------------------  INPUTS TO MAIN PROGRAM  -------------------
#
# DogboneAnalytic_ML_subcycle, cut short and asserting the coarse-fine metric identity rather
# than any answer. Advection_ML_cf_metrics checks the same identity at a refinement ratio of
# 2, where three fine faces never have to tile one coarse face; this is the ratio-3 grid the
# DogboneAnalytic_ML golds and the volume-conservation lanes are all measured on, so it is the
# one whose tiling those numbers actually rest on. check_cf_metrics aborts past
# remora.check_cf_tol, so completing the run is the assertion.
#
remora.prob_name = DogboneAnalytic

remora.max_step = 4

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
remora.do_substep = 1

amr.max_level       = 1       # maximum level number allowed
amr.ref_ratio_vect = 3 3  1

# COARSE-FINE METRIC CHECK
# Assert that the fine cell edges sum to the coarse edge across the interface. Everything
# set_2d_cf_bcs does rests on that identity: it hands the child the parent's mass flux per unit
# edge length and lets each fine face multiply its own length back in. At a ratio of 3 the sum
# is over three faces rather than two, so it is the case where an off-by-one in the stencil or
# a metric rescaled instead of summed would show.
remora.check_cf_metrics = 1
# The default tolerance, 1e-12, is four orders above what any of these grids measures
# (0, 1.4e-16, 2.5e-16). 1e-14 keeps the assertion within a factor of 40 of the answer.
remora.check_cf_tol     = 1e-14

# DIAGNOSTICS & VERBOSITY
remora.sum_interval  = 1       # timesteps between integrated/max quantities, if remora.v > 0
remora.v             = 0       # verbosity in REMORA.cpp (0: none, 1: integrated quantities, etc, 2: print boxes)

# CHECKPOINT FILES
remora.check_file      = chk        # root name of checkpoint file
remora.check_int       = -57600     # number of timesteps between checkpoints

# PLOTFILES
remora.plot_file     = plt_ml        # prefix of plotfile name
remora.plot_int      = -1            # nothing compares plotfiles here
remora.plot_vars_3d  = salt temp x_velocity y_velocity z_velocity
remora.plotfile_type = amrex
remora.expand_plotvars_to_unif_rr = 1

# SOLVER CHOICE
remora.tracer_horizontal_advection_scheme = "upstream3" # upstream3 or centered4

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

# PROBLEM PARAMETERS (optional)
remora.R0    = 1027.0  # background density value (Kg/m3) used in Linear Equation of State
remora.S0    = 35.0    # background salinity (nondimensional) constant
remora.T0    = 10.0    # background potential temperature (Celsius) constant
remora.Tcoef = 1.7e-4  # linear equation of state parameter (1/Celsius)
remora.Scoef = 7.6e-4     # linear equation of state parameter (nondimensional)
remora.rho0  = 1025.0  # Mean density (Kg/m3) used when Boussinesq approx is inferred

# These files can be found at https://github.com/seahorce-scidac/REMORA-data
remora.ic_type       = "analytic"

remora.refinement_indicators = velup veldown
remora.coupling_type = "TwoWay"
remora.velup.max_level=1
remora.velup.field_name = x_velocity
remora.velup.value_greater = 0.05
remora.veldown.max_level=1
remora.veldown.field_name = x_velocity
remora.veldown.value_less = -0.05
