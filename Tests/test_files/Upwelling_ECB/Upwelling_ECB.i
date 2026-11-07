# Synthetic integration test; parameters exercise clean ROMS ECB, not a calibrated estuary.
# ------------------  INPUTS TO MAIN PROGRAM  -------------------
remora.prob_name = Upwelling

remora.max_step = 3

remora.omp_tile_size = 1024 1024 1024

# PROBLEM SIZE & GEOMETRY
remora.prob_lo     =      0.     0.    -150.
remora.prob_hi     =  41000. 80000.       0.

remora.n_cell           =  8      8       8

remora.is_periodic = 1 0 0

remora.bc.ylo.type = "SlipWall"
remora.bc.yhi.type = "SlipWall"

# TIME STEP CONTROL
remora.fixed_dt       = 300.0 # Timestep size (seconds)

remora.ndtfast  = 30 # Ratio of baroclinic to barotropic time step

# DIAGNOSTICS & VERBOSITY
remora.sum_interval  = 1       # timesteps between integrated/max quantities, if remora.v > 0
remora.v             = 0       # verbosity in REMORA.cpp (0: none, 1: integrated quantities, etc, 2: print boxes)

# CHECKPOINT FILES
remora.check_file      = chk        # root name of checkpoint file
remora.check_int       = 2      # number of timesteps between checkpoints

# PLOTFILES
remora.plot_file     = plt        # prefix of plotfile name
remora.plot_int      = 1        # number of timesteps between plotfiles
remora.plot_vars_3d  = salt temp tracer x_velocity y_velocity z_velocity biology
remora.plotfile_type = amrex

# SOLVER CHOICE
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


# ECB with explicit reference-unit parameters.
remora.biology_model = ecb
remora.biology_ic_type = constant
remora.nscalar = 1
remora.bulk_fluxes = true
remora.surface_radiation_flux = 150.0
remora.uwind = 5.0
remora.vwind = 0.0
remora.ecb.BioIter = 2
remora.ecb.AttSW = 0.04
remora.ecb.AttChl = 0.02486
remora.ecb.Chl2C_m = 0.0535
remora.ecb.ChlMin = 0.001
remora.ecb.CoagR = 0.015
remora.ecb.D_p5NH4 = 0.1
remora.ecb.I_thNH4 = 0.0095
remora.ecb.K_NH4 = 6.7
remora.ecb.K_NO3 = 2.0
remora.ecb.K_PO4 = 32.0
remora.ecb.K_Phy = 2.0
remora.ecb.LDeRRN = 0.09
remora.ecb.LDeRRC = 0.09
remora.ecb.NitriR = 0.05
remora.ecb.PARfrac = 0.43
remora.ecb.PhyCN = 6.625
remora.ecb.R_P2N = 0.0625
remora.ecb.PhyIS = 0.08
remora.ecb.PhyMin = 0.001
remora.ecb.PhyMR = 0.05
remora.ecb.SDeRRN = 0.09
remora.ecb.SDeRRC = 0.09
remora.ecb.RDeRRN = 0.01
remora.ecb.RDeRRC = 0.01
remora.ecb.Vp0 = 0.55
remora.ecb.wLDet = 5.0
remora.ecb.wPhy = 0.1
remora.ecb.wSDet = 0.1
remora.ecb.ZooAE_N = 0.75
remora.ecb.ZooBM = 0.1
remora.ecb.ZooCN = 6.625
remora.ecb.ZooER = 0.1
remora.ecb.ZooGR = 0.03
remora.ecb.ZooMin = 0.001
remora.ecb.ZooMR = 0.025
remora.ecb.ElDON = 0.2
remora.ecb.EsDON = 0.1
remora.ecb.deltN = 0.15
remora.ecb.deltC = 0.275
remora.ecb.gammaC = 0.15
remora.ecb.slCexc = 0.45
remora.ecb.pCO2air = 400.0
# Constant initial concentrations, set by tracer name.
remora.biology_initial_NO3           = 8.0
remora.biology_initial_NH4           = 0.5
remora.biology_initial_chlorophyll   = 1.0
remora.biology_initial_phytoplankton = 0.3
remora.biology_initial_zooplankton   = 0.1
remora.biology_initial_LdetritusN    = 0.2
remora.biology_initial_SdetritusN    = 0.1
remora.biology_initial_DON           = 3.0
remora.biology_initial_PO4           = 0.5
remora.biology_initial_LdetritusC    = 1.2
remora.biology_initial_SdetritusC    = 0.6
remora.biology_initial_TIC           = 2100.0
remora.biology_initial_alkalinity    = 2300.0
remora.biology_initial_DOC           = 50.0
remora.biology_initial_oxygen        = 200.0
remora.biology_initial_Rdon          = 1.0
remora.biology_initial_Rdoc          = 40.0
remora.biology_initial_calcium       = 10000.0

# Fill new arrays with signaling NaN and trap invalid operations, so a biology
# kernel that reads scratch it has not written fails instead of using garbage.
amrex.init_snan = 1
amrex.fpe_trap_invalid = 1
