# ------------------  INPUTS TO MAIN PROGRAM  -------------------
max_step = 10
stop_time = 30000.0

amrex.fpe_trap_invalid = 0

fabarray.mfiter_tile_size = 16 16 1024

# PROBLEM SIZE & GEOMETRY
geometry.prob_lo     =      0.     0.    -150.
geometry.prob_hi     =  41000. 41000.       0.

amr.n_cell           = 81     81      16

# periodic in x to match WRF setup
geometry.is_periodic = 1 1 0
#ylo.type = "SlipWall"
#yhi.type = "SlipWall"
zlo.type = "SlipWall"
zhi.type = "SlipWall"

# TIME STEP CONTROL
romsx.fixed_dt       = 300.0 # Timestep size (seconds)
# NDTFAST  = 30.0 # Number of baratropic steps => 300.0/30.0 = 10.0
romsx.fixed_fast_dt  = 10.0 # Baratropic timestep size (seconds)
# romsx.fixed_fast_dt  = 300.0 # Baratropic timestep size (seconds) testing value
romsx.fixed_ndtfast_ratio  = 30 # Baratropic timestep size (seconds)

# DIAGNOSTICS & VERBOSITY
romsx.sum_interval   = 10       # timesteps between computing mass
romsx.v              = 0       # verbosity in ROMSX.cpp (0: none, 1: print boxes, etc, 2: print values)
amr.v                = 1       # verbosity in Amr.cpp

# REFINEMENT / REGRIDDING
amr.max_level       = 0       # maximum level number allowed

# CHECKPOINT FILES
romsx.check_file      = chk        # root name of checkpoint file
romsx.check_int       = -57600      # number of timesteps between checkpoints

# PLOTFILES
romsx.plot_file_1     = plt        # prefix of plotfile name
romsx.plot_int_1      = 100          # number of timesteps between plotfiles
romsx.plot_vars_1     = salt temp scalar x_velocity y_velocity z_velocity
romsx.plotfile_type   = amrex

# SOLVER CHOICE
romsx.use_coriolis = false
romsx.flat_bathymetry=true
romsx.horizontal_advection_scheme = "upstream3" # upstream3 or centered4
romsx.spatial_order = 2

# Linear EOS parameters
romsx.R0    = 1027.0  # background density value (Kg/m3) used in Linear Equation of State
romsx.S0    = 35.0    # background salinity (nondimensional) constant
romsx.T0    = 14.0    # background potential temperature (Celsius) constant
romsx.Tcoef = 1.7e-4  # linear equation of state parameter (1/Celcius)
romsx.Scoef = 0.0     # linear equation of state parameter (nondimensional)
romsx.rho0  = 1025.0  # Mean density (Kg/m3) used when Boussinesq approx is inferred

# Coriolis params
romsx.coriolis_f0 = 0.0
romsx.coriolis_beta = 0.0

# PROBLEM PARAMETERS (velocity fields)
prob.u_0   =1.0e-1
prob.v_0   = -1.0e-1
