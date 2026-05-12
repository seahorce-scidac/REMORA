# ------------------  INPUTS TO MAIN PROGRAM  -------------------
remora.max_step = 0
remora.stop_time = 30000.0

amrex.fpe_trap_invalid = 1

# PROBLEM SIZE & GEOMETRY
remora.prob_lo     =      0.     0.    -150.
remora.prob_hi     =  41000. 80000.       0.

remora.n_cell           =  41     80      16

remora.is_periodic = 1 1 0

# TIME STEP CONTROL
remora.fixed_dt       = 300.0
remora.fixed_fast_dt  = 10.0
remora.fixed_ndtfast_ratio  = 30

# DIAGNOSTICS & VERBOSITY
remora.sum_interval   = 1
remora.v              = 0
amr.v                 = 1

# CHECKPOINT FILES
remora.check_file      = chk
remora.check_int       = -57600

# PLOTFILES
remora.plot_file     = plt
remora.plot_int      = 100
remora.plot_vars_3d  = salt temp x_velocity y_velocity z_velocity
remora.plotfile_type = amrex

# SOLVER CHOICE
remora.use_coriolis = true
remora.flat_bathymetry=true
remora.tracer_horizontal_advection_scheme = "upstream3"
remora.spatial_order = 2

# Linear EOS parameters
remora.R0    = 1027.0
remora.S0    = 35.0
remora.T0    = 14.0
remora.Tcoef = 1.7e-4
remora.Scoef = 0.0
remora.rho0  = 1025.0

# Coriolis params
remora.coriolis_f0 = -8.26e-5
remora.coriolis_beta = 0.0

# PROBLEM PARAMETERS (velocity fields)
remora.prob.u_0   = 0.0
remora.prob.v_0   = 0.0
remora.prob.z0    = 0.1
remora.prob.zRef  = 80.0e-3
remora.prob.uRef  = 8.0e-3
