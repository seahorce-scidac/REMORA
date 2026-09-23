# ------------------  INPUTS TO MAIN PROGRAM  -------------------
# Channel_Test with a refined patch against the periodic x boundary.
#
# The patch covers x = 0..50000 of a domain periodic in x, so its low edge sits on the seam
# and its high edge is an ordinary coarse-fine interface in open water. That geometry is what
# exercises the periodic-boundary handling in BuildMask, fill_ghost_kcomps, FillCoarsePatchMap
# and the ubar/vbar average-down: with any of them reverted the step-1 answer changes.
# Advection_ML's patch never touches its periodic edges, so this lane is the only one that does.
#
# It is also the only multi-level lane with Coriolis on, and its gold comes from one rank with
# one box per level. Under MPI the level-1 box splits at j = 30, and the first version of this
# lane found that the answer depended on that split: set_2d_cf_bcs wrote the interface faces
# but not the neighbouring box's ghost copies, and Coriolis read them. Both ML Dogbone inputs
# run Coriolis off, which is why they never saw it.
remora.prob_name = ChannelTest

remora.max_step = 10

amrex.fpe_trap_invalid = 1

# PROBLEM SIZE & GEOMETRY
remora.prob_lo     =      0.      0.     -50.
remora.prob_hi     = 100000. 300000.       0.

remora.n_cell           =  20     60      50

remora.is_periodic = 1 0 0

remora.bc.ylo.type = "SlipWall"
remora.bc.yhi.type = "SlipWall"

# TIME STEP CONTROL
remora.fixed_dt       = 400.0 # Timestep size (seconds)

remora.ndtfast = 10

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
remora.vertical_mixing_type = gls

remora.gls_P = 3.0
remora.gls_M = 1.5
remora.gls_N = -1.0
remora.gls_Kmin = 7.6e-6
remora.gls_Pmin = 1.0e-12

remora.gls_cmu0 = 0.5477
remora.gls_c1 = 1.44
remora.gls_c2 = 1.92
remora.gls_c3m = -0.4
remora.gls_c3p = 1.0
remora.gls_sigk = 1.0
remora.gls_sigp = 1.3

remora.Zob = 0.002
remora.Zos = 0.002

# turbulence closure parameters
remora.Akk_bak = 5.0e-6
remora.Akp_bak = 5.0e-6
remora.Akv_bak = 1.0e-5

# Linear EOS parameters
remora.R0    = 1027.0  # background density value (Kg/m3) used in Linear Equation of State
remora.S0    = 15.0    # background salinity (nondimensional) constant
remora.T0    = 10.0    # background potential temperature (Celsius) constant
remora.Tcoef = 1.7e-4  # linear equation of state parameter (1/Celsius)
remora.Scoef = 7.6e-4     # linear equation of state parameter (nondimensional)
remora.rho0  = 1025.0  # Mean density (Kg/m3) used when Boussinesq approx is inferred

remora.tcline = 25.0

# Coriolis params
remora.use_coriolis = true
remora.coriolis_type = beta_plane
remora.coriolis_f0 = 1.0e-4
remora.coriolis_beta = 0.0

# REFINEMENT
# Half the periodic width from x = 0, and half the channel from y = 0. The low x edge is the
# seam; the high x edge and the y = 150000 edge are coarse-fine interfaces.
remora.refinement_indicators = seam
remora.coupling_type = "TwoWay"
remora.seam.max_level = 1
remora.seam.in_box_lo =     0.      0. -50.
remora.seam.in_box_hi = 50000. 150000.   0.

amr.max_level      = 1
amr.ref_ratio_vect = 2 2 1
