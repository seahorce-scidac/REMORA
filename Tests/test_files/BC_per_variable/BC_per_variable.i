# ------------------  INPUTS TO MAIN PROGRAM  -------------------
#
# Per-variable boundary specification, pinning the West, South, East, North order: all four
# sides carry a different condition. BC_per_side.i sets the same four faces by side, and the
# two runs must agree -- see add_test_equiv in Tests/CTestList.cmake.
#
# The dye is zero in the initial condition and enters only through the western inflow value,
# so a nonzero tracer field means that value was read.
#
remora.prob_name = DoubleGyre

remora.max_step = 10

amrex.fpe_trap_invalid = 1

# PROBLEM SIZE & GEOMETRY
remora.prob_lo     =      0.      0.     -500.
remora.prob_hi     = 1000000. 2000000.       0.

remora.n_cell           =  54     108      4

remora.is_periodic = 0 0 0

remora.nscalar = 1

# BOUNDARY CONDITIONS
#                        West      South       East      North
remora.bc.temp.type   =  inflow    noslipwall  outflow   slipwall
remora.bc.salt.type   =  inflow    noslipwall  outflow   slipwall
remora.bc.tracer.type =  inflow    noslipwall  outflow   slipwall
remora.bc.u.type      =  inflow    noslipwall  outflow   slipwall
remora.bc.v.type      =  inflow    noslipwall  outflow   slipwall
remora.bc.w.type      =  inflow    noslipwall  outflow   slipwall
remora.bc.ubar.type   =  inflow    noslipwall  outflow   slipwall
remora.bc.vbar.type   =  inflow    noslipwall  outflow   slipwall
remora.bc.zeta.type   =  inflow    noslipwall  outflow   slipwall
remora.bc.tke.type    =  inflow    noslipwall  outflow   slipwall

remora.boundary_per_variable = 1

# One entry per variable, applying to every side of that variable that takes a Dirichlet value --
# here the western inflow and the southern noslipwall. The velocity is nonzero so that the
# noslipwall face is pinned to a real tangential value rather than agreeing with the other lane
# only because both are zero; the normal component is zeroed there by init_bcs in both lanes, and
# BC_per_side.i states the same entry on ylo so the two stay equivalent.
#
# The inflow is diffusive, not advective: ubar and vbar take no inflow value (uses_velocity_input
# covers the 3D velocities alone), so the barotropic normal velocity at the face stays zero and
# the correction holds the net flux near zero. Raising this entry tenfold moves the dye ~7%.
remora.bc.temp.value    = 8.0
remora.bc.salt.value    = 35.5
remora.bc.tracer.value  = 1.0
remora.bc.u.velocity    = 0.1 0. 0.
remora.bc.v.velocity    = 0.1 0. 0.
remora.bc.w.velocity    = 0.1 0. 0.

# TIME STEP CONTROL
remora.fixed_dt       = 3600.0 # Timestep size (seconds)

remora.fixed_ndtfast_ratio = 20

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
remora.plot_vars_3d  = salt temp tracer x_velocity y_velocity z_velocity
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
remora.tnu2_scalar = 1280.0

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
