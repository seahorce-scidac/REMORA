# ------------------  SeaIceBox  ------------------
# A small closed basin with the myROMS sea-ice model: water just above freezing
# under a cold, windy, snowing atmosphere, with an ice pack over its western
# 60% (analytic initial state). Exercises the ice bulk fluxes, Mellor-Kantha
# thermodynamics (including snow), EVP momentum and MPDATA advection. The
# ice drives the ocean through the surface stress and heat and salt fluxes,
# which is what the plotfile comparison sees.

remora.prob_name = SeaIceBox

max_step = 10

remora.prob_lo     =      0.      0.   -50.
remora.prob_hi     =  64000.  32000.     0.
remora.n_cell      =     32      16      8

remora.is_periodic = 0 0 0
remora.bc.xlo.type = "SlipWall"
remora.bc.xhi.type = "SlipWall"
remora.bc.ylo.type = "SlipWall"
remora.bc.yhi.type = "SlipWall"

remora.fixed_dt = 300.0
remora.ndtfast  = 30

remora.v = 0
remora.sum_interval = 1

remora.check_file = chk
remora.check_int  = -1

remora.plot_file     = plt
remora.plot_int      = 10
remora.plotfile_type = amrex
remora.plot_vars_3d  = salt temp x_velocity y_velocity z_velocity
remora.plot_vars_2d  = seaice zeta

remora.theta_s = 3.0
remora.theta_b = 0.4
remora.tcline  = 20.0

remora.eos_type = nonlinear
remora.rho0     = 1025.0
remora.g        = 9.81

remora.use_coriolis  = true
remora.coriolis_type = beta_plane
remora.coriolis_f0   = 1.4e-4
remora.coriolis_beta = 0.0

remora.horizontal_mixing_type = constant
remora.visc2     = 5.0
remora.tnu2_temp = 5.0
remora.tnu2_salt = 5.0

remora.bottom_stress_type = quadratic
remora.rdrag2 = 3.0e-3

remora.vertical_mixing_type = GLS
remora.Akv_bak = 1.0e-5
remora.Akt_bak = 1.0e-6

# Constant atmosphere; negative rain is snow
remora.bulk_fluxes = true
remora.bulk_flux_scheme = coare35
remora.eminusp     = true
remora.uwind_type  = constant
remora.vwind_type  = constant
remora.uwind = 6.0
remora.vwind = -4.0
remora.air_temperature = -15.0
remora.air_humidity    = 0.7
remora.air_pressure    = 1013.48
remora.surface_radiation_flux = 50.0
remora.cloud = 0.5
remora.rain  = -2.0e-5

remora.ice_model = mellor_kantha
ice.mom_bulk = true
ice.convsnow = true
ice.nEVP     = 30
ice.Pstar    = 25000.0
ice.bc_type  = closed closed closed closed
