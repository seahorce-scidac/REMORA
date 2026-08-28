# ------------------  INPUTS TO MAIN PROGRAM  -------------------
# Upwelling carrying the Fennel biology model. Physics and biology are both analytic, so
# this case needs no NetCDF and runs anywhere the other regression cases do.
#
# The option set below is deliberately wide: carbon, oxygen, the non-conservative
# alkalinity budget, the river DON pool, the dated atmospheric pCO2 and the Wanninkhof
# (2014) gas transfer are all on at once, so one gold file covers all of them. See the
# BIOLOGY block for what each option pulls in and why these values were chosen.
#
# remora.nscalar = 1 is deliberate. Upwelling initializes its dye to zero everywhere, so
# the tracer column is not testing dye physics -- it pins the component layout. Biology
# starts at Bio_comp = Tracer_comp + nscalar, so if that offset were wrong a biology
# tracer would land in the dye slot and the tracer column would stop being zero.
remora.prob_name = Upwelling

remora.max_step = 10

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

remora.ndtfast  = 30 # Ratio of baroclinic to barotropic time step

# DIAGNOSTICS & VERBOSITY
remora.sum_interval  = 1       # timesteps between computing mass
remora.v             = 0       # verbosity in REMORA.cpp (0: none, 1: print boxes, etc, 2: print values)
amr.v                = 1       # verbosity in Amr.cpp

# CHECKPOINT FILES
remora.check_file      = chk        # root name of checkpoint file
remora.check_int       = -57600      # number of timesteps between checkpoints

# PLOTFILES
remora.plot_file     = plt        # prefix of plotfile name
remora.plot_int      = 100        # number of timesteps between plotfiles
remora.plot_vars_3d  = salt temp tracer x_velocity y_velocity z_velocity fennel
remora.plotfile_type = amrex

# BIOLOGY
remora.biology_model    = fennel
remora.nscalar          = 1
remora.biology_ic_type  = analytic

# Carbon and oxygen carry the tracers the remaining options act on: alkalinity and TIC
# for the two CO2 options, dissolved oxygen for the O2 Schmidt number. Neither is a
# prerequisite in the input parser, but with both off the four options below would have
# nothing to change and the test would pass no matter what they did.
remora.fennel.carbon    = true
remora.fennel.oxygen    = true

# Alkalinity responds to nitrate uptake, nitrification, zooplankton metabolism and
# remineralization instead of being diagnosed from salinity (ROMS TALK_NONCONSERV).
remora.fennel.talk_nonconserv = true

# River dissolved organic nitrogen, plus its carbon counterpart since carbon is on
# (ROMS RIVER_DON). The analytic initial condition seeds both pools; a zero start would
# leave the option a no-op. Neither pool sinks, matching ROMS.
remora.fennel.river_don = true
remora.fennel.RDeRRN    = 0.03
remora.fennel.RDeRRC    = 0.03

# Atmospheric pCO2 from the secular trend rather than a constant (ROMS PCO2AIR_SECULAR).
# The trend is anchored to 1951 and goes negative long before it, so the run needs a real
# date: remora.time_ref below places model time zero at 1 March 2020. That date is chosen
# to fall after 29 February of a leap year, so the day-of-year the calendar hands to the
# trend depends on the leap-day branch of REMORA_DateClock.H.
remora.fennel.pco2air_type = secular

# Wanninkhof (2014) Schmidt numbers and rate coefficient for both gases, in place of the
# Wanninkhof (1992) default (ROMS RW14_CO2_SC and RW14_OXYGEN_SC). The OCMIP oxygen set
# is the one formulation of the three that no regression test exercises.
remora.fennel.co2_schmidt    = rw14
remora.fennel.oxygen_schmidt = rw14

# Denitrification stays off, which is not just the default: the alkalinity increment from
# the sediment remineralization return lives in the branch ROMS takes when BIO_SEDIMENT is
# on and DENITRIFICATION is off, so turning denitrification on here would leave that
# increment untested.
remora.fennel.denitrification = false

# Reference date and calendar for the model clock, dating the pCO2 trend above.
remora.time_ref = 20200301.0

# Fennel needs shortwave radiation, which is only allocated under bulk fluxes.
remora.bulk_fluxes             = true
remora.surface_radiation_flux  = 300.0   # W/m2
remora.uwind                   = 5.0
remora.vwind                   = 0.0

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
