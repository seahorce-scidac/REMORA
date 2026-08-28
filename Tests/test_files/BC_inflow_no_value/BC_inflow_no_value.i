# ------------------  INPUTS TO MAIN PROGRAM  -------------------
#
# An inflow face with no inflow value for temperature. Nothing else supplies one, so the run
# must stop at setup naming the missing input, rather than write the placeholder default into
# the ghost cells. Velocity inflow has always insisted on a value the same way.
#
remora.prob_name = DoubleGyre

remora.max_step = 10

amrex.fpe_trap_invalid = 1

# PROBLEM SIZE & GEOMETRY
remora.prob_lo     =      0.      0.     -500.
remora.prob_hi     = 1000000. 2000000.       0.

remora.n_cell           =  54     108      4

remora.is_periodic = 0 0 0

# BOUNDARY CONDITIONS
remora.bc.xlo.type = inflow       # West
remora.bc.ylo.type = noslipwall   # South
remora.bc.xhi.type = outflow      # East
remora.bc.yhi.type = slipwall     # North

remora.boundary_per_variable = 0

# Velocity and salt are given; temp is not, which is the point.
remora.bc.xlo.salt     = 35.5
remora.bc.xlo.velocity = 0. 0. 0.

# TIME STEP CONTROL
remora.fixed_dt       = 3600.0 # Timestep size (seconds)

remora.ndtfast = 20

# DIAGNOSTICS & VERBOSITY
remora.sum_interval  = 1       # timesteps between integrated/max quantities, if remora.v > 0
remora.v             = 0       # verbosity in REMORA.cpp (0: none, 1: integrated quantities, etc, 2: print boxes)
amr.v                = 1       # verbosity in Amr.cpp

# CHECKPOINT FILES
remora.check_file      = chk        # root name of checkpoint file
remora.check_int       = -57600      # number of timesteps between checkpoints

# PLOTFILES
remora.plot_file     = plt        # prefix of plotfile name
remora.plot_int      = 10         # number of timesteps between plotfiles
remora.plot_vars_3d  = salt temp x_velocity y_velocity z_velocity
remora.plotfile_type = amrex

# SOLVER CHOICE
remora.tracer_horizontal_advection_scheme = "upstream3" # upstream3 or centered4

remora.Zob = 0.02
remora.Zos = 0.02

remora.rdrag = 8.0e-7

# turbulence closure parameters
remora.Akk_bak = 5.0e-6
remora.Akp_bak = 5.0e-6
remora.Akv_bak = 1.0
remora.Akt_bak = 1.0

remora.theta_s = 0.0
remora.theta_b = 0.0
remora.tcline = 1e16

remora.horizontal_mixing_type = constant
remora.visc2 = 1280.0
remora.tnu2_salt = 1280.0
remora.tnu2_temp = 1280.0

# Linear EOS parameters
remora.R0    = 1028.0  # background density value (Kg/m3) used in Linear Equation of State
remora.S0    = 35.0    # background salinity (nondimensional) constant
remora.T0    = 5.0    # background potential temperature (Celsius) constant
remora.Tcoef = 1.0e-4  # linear equation of state parameter (1/Celsius)
remora.Scoef = 7.6e-4     # linear equation of state parameter (nondimensional)
remora.rho0  = 1025.0  # Mean density (Kg/m3) used when Boussinesq approx is inferred

# Coriolis params
remora.use_coriolis = true
remora.coriolis_type = beta_plane
remora.coriolis_f0 = 7.3e-5
remora.coriolis_beta = 2.0e-11
