# ------------------  INPUTS TO MAIN PROGRAM  -------------------
#
# The subcycled driver at a timestep ratio of 1 must reproduce the lockstep driver.
#
# Run twice from CTestList: remora.do_substep=0 goes through timeStepML, and
# remora.do_substep=1 remora.dt_ref_ratio=1 goes through the recursive timeStep with every level
# on dt[0]. The two plotfiles must agree, which separates "the recursive driver changed
# something" from "subcycling changes the answer" -- the second is expected, the first is a
# bug. No gold file: neither run is the baseline.
#
# fixed_dt is half Advection_ML's. At dt = 100 with dt_ref_ratio = 1, level 1 runs the
# coarse barotropic step on a twice-finer grid and goes unstable, which would test the CFL
# limit rather than the drivers.
#
# This holds with set_2d_cf_bcs active only because Advection has a flat bottom: D is
# the same either side of the interface, so the mass-flux form reduces to the interpolation it
# replaces. On varying bathymetry the two differ and ratio 1 would no longer match lockstep --
# correctly, since the 2D coupling is then genuinely different. See DogboneAnalytic_ML_subcycle.
#
# Refluxing is not degenerate here and is switched off in the comparison: it is a correction
# the lockstep driver never applies, so leaving it on would test a feature rather than the
# drivers. It moves the tracer by 3e-4, and only the tracer -- temp and salt are uniform, so
# their flux mismatch across the interface is zero.
#
remora.prob_name = Advection

remora.max_step = 20
remora.stop_time = 300000.0

amrex.fpe_trap_invalid = 1

# PROBLEM SIZE & GEOMETRY
remora.prob_lo     =      0.     0.    -150.
remora.prob_hi     =  40000. 40000.       0.

remora.n_cell           = 80     80      16

remora.is_periodic = 1 1 0

# TIME STEP CONTROL
remora.fixed_dt       = 50.0 # Timestep size (seconds)
remora.ndtfast  = 10

# DIAGNOSTICS & VERBOSITY
remora.sum_interval   = 20       # timesteps between computing mass
remora.v              = 0       # verbosity in REMORA.cpp (0: none, 1: print boxes, etc, 2: print values)
amr.v                = 1       # verbosity in Amr.cpp

# REFINEMENT / REGRIDDING
amr.max_level       = 1       # maximum level number allowed
amr.ref_ratio_vect = 2 2 1
amr.regrid_int=1

# CHECKPOINT FILES
remora.check_file      = chk        # root name of checkpoint file
remora.check_int       = -57600      # number of timesteps between checkpoints

# PLOTFILES
remora.plot_file     = plt        # prefix of plotfile name
remora.plot_int      = 20         # number of timesteps between plotfiles
# Dye is opt-in (remora.nscalar defaults to 0), and this case advects a dye blob, so ask for one.
remora.nscalar        = 1
remora.plot_vars_3d  = salt temp tracer x_velocity y_velocity z_velocity
remora.plotfile_type = amrex
remora.expand_plotvars_to_unif_rr = 1

# SOLVER CHOICE
remora.use_coriolis = false
remora.tracer_horizontal_advection_scheme = "centered4" # upstream3 or centered4

# Linear EOS parameters
remora.R0    = 1027.0  # background density value (Kg/m3) used in Linear Equation of State
remora.S0    = 35.0    # background salinity (nondimensional) constant
remora.T0    = 14.0    # background potential temperature (Celsius) constant
remora.Tcoef = 0.0 #1.7e-4  # linear equation of state parameter (1/Celsius)
remora.Scoef = 0.0 #1.0e-4  # linear equation of state parameter (nondimensional)
remora.rho0  = 1025.0  # Mean density (Kg/m3) used when Boussinesq approx is inferred

# Coriolis params
remora.coriolis_f0 = 0.0
remora.coriolis_beta = 0.0

remora.rdrag=0.0

# PROBLEM PARAMETERS (velocity fields)
remora.prob.u_0   =1.0e-0
remora.prob.v_0   = -1.0e-0

remora.refinement_indicators = scalar
remora.scalar.max_level = 1
remora.scalar.field_name = tracer
remora.scalar.value_greater = 0.5
remora.scalar.start_time = 200
remora.coupling_type = "TwoWay"
