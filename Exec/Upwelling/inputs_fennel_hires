# ------------------  INPUTS TO MAIN PROGRAM  -------------------
# Upwelling carrying the Fennel biology model. Physics and biology are both analytic, so
# this case needs no NetCDF and runs anywhere the other regression cases do.
#
# remora.nscalar = 1 is deliberate. Upwelling initializes its dye to zero everywhere, so
# the tracer column is not testing dye physics -- it pins the component layout. Biology
# starts at Bio_comp = Tracer_comp + nscalar, so if that offset were wrong a biology
# tracer would land in the dye slot and the tracer column would stop being zero.
remora.prob_name = Upwelling

remora.max_step = 0

amrex.fpe_trap_invalid = 1

# PROBLEM SIZE & GEOMETRY
remora.prob_lo     =      0.     0.    -150.
remora.prob_hi     =  41000. 80000.       0.

remora.n_cell           =  41     80      16

remora.is_periodic = 1 0 0

remora.bc.ylo.type = "SlipWall"
remora.bc.yhi.type = "SlipWall"

# TIME STEP CONTROL
remora.fixed_dt       = 300.0 # Timestep size (seconds)

remora.fixed_fast_dt  = 10.0 # Baratropic timestep size (seconds)

remora.fixed_ndtfast_ratio  = 30 # Ratio of baroclinic to barotropic time step

# DIAGNOSTICS & VERBOSITY
remora.sum_interval  = 1       # timesteps between computing mass
remora.v             = 0       # verbosity in REMORA.cpp (0: none, 1: print boxes, etc, 2: print values)
amr.v                = 1       # verbosity in Amr.cpp

# CHECKPOINT FILES
remora.check_file      = chk        # root name of checkpoint file
remora.check_int       = -57600      # number of timesteps between checkpoints

# PLOTFILES
remora.plot_file     = plt        # prefix of plotfile name
remora.plot_int      = 1          # number of timesteps between plotfiles
remora.plot_vars_3d  = salt temp tracer x_velocity y_velocity z_velocity fennel
remora.plotfile_type = amrex

# BIOLOGY
remora.biology_model    = fennel
remora.nscalar          = 1
remora.biology_ic_type  = analytic

# Fennel needs shortwave radiation, which is only allocated under bulk fluxes.
remora.bulk_fluxes             = true
remora.surface_radiation_flux  = 300.0   # W/m2
remora.uwind                   = 5.0
remora.vwind                   = 0.0

# SOLVER CHOICE
remora.flat_bathymetry = false
remora.tracer_horizontal_advection_scheme = "upstream3" # upstream3 or centered4

# Linear EOS parameters
remora.R0    = 1027.0  # background density value (Kg/m3) used in Linear Equation of State
remora.S0    = 35.0    # background salinity (nondimensional) constant
remora.T0    = 14.0    # background potential temperature (Celsius) constant
remora.Tcoef = 1.7e-4  # linear equation of state parameter (1/Celsius)
remora.Scoef = 0.0     # linear equation of state parameter (nondimensional)
remora.rho0  = 1025.0  # Mean density (Kg/m3) used when Boussinesq approx is inferred

remora.tcline = 25.0

# Coriolis params
remora.use_coriolis = true
remora.coriolis_type = beta_plane
remora.coriolis_f0 = -8.26e-5
remora.coriolis_beta = 0.0

# HIGH-RESOLUTION BATHYMETRY WITH BIOLOGY -- exact value assertions, no baseline
# Bathymetry comes from level 1 and is averaged down while the Fennel tracers are initialized
# by the analytic profile. Six of the seven Fennel tracers are constants in
# REMORA_InitAnalyticBiology_BioToy.H, so their expected level-0 values are known outright
# rather than snapshotted -- and because the dye sits at Tracer_comp with biology starting at
# Bio_comp = Tracer_comp + nscalar, asserting tracer == 0 alongside them is what catches a
# component-offset regression that would shift a biology tracer into the dye slot.
# max_step = 0 makes this an initialize-and-dump run: ~1 s, and no time stepper in the way.
amr.max_level                     = 1
amr.ref_ratio_vect                = 3 3 1
remora.expand_plotvars_to_unif_rr = 1
remora.hires_grid_level           = 1
