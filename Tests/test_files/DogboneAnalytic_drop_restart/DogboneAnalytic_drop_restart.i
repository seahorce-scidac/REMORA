# ------------------  INPUTS TO MAIN PROGRAM  -------------------
# A run with dropped level-0 boxes, restarted from its step-5 checkpoint, must finish
# where the continuous run does. The restart keeps the checkpoint's BoxArray.
remora.prob_name = DogboneAnalytic

remora.max_step = 10

# PROBLEM SIZE & GEOMETRY
remora.prob_lo     =   0.       0.  -10.
remora.prob_hi     =   8400.  750.       0.

remora.n_cell           =  128 64 8

amr.max_grid_size_x   = 16
amr.max_grid_size_y   = 16
amr.blocking_factor_x = 16
amr.blocking_factor_y = 16

remora.is_periodic = 0 0 0

remora.bc.xlo.type = "slipwall"
remora.bc.xhi.type = "slipwall"
remora.bc.ylo.type = "slipwall"
remora.bc.yhi.type = "slipwall"

# TIME STEP CONTROL
remora.fixed_dt            = 6.0 # Timestep size (seconds)
remora.ndtfast = 20

# DIAGNOSTICS & VERBOSITY
remora.sum_interval  = -1
remora.v             = 0

# CHECKPOINT FILES
remora.check_file      = chk
remora.check_int       = 5

# PLOTFILES
remora.plot_file     = plt
remora.plot_int      = 10
remora.plot_vars_3d  = salt temp x_velocity y_velocity z_velocity
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

# PROBLEM PARAMETERS
remora.R0    = 1027.0
remora.S0    = 35.0
remora.T0    = 10.0
remora.Tcoef = 1.7e-4
remora.Scoef = 7.6e-4
remora.rho0  = 1025.0

remora.ic_type       = "analytic"

# Land from x = 1400 to 7000 outside the channel, and on the boundary ring
remora.prob.mask_x_lo = 1400.
remora.prob.mask_x_hi = 7000.
remora.prob.mask_ring = true

# Horizontal mixing reaches across box edges into the land
remora.horizontal_mixing_type = constant
remora.harmonic_mixing_type   = geopotential
remora.visc2                  = 5.0
remora.tnu2_temp              = 5.0
remora.tnu2_salt              = 5.0
remora.drop_dry_boxes = true
