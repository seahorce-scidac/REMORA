# A zonal wind stress over a periodic, flat-bottom ocean of depth 500 m.
# With no drag, viscosity or Coriolis, u = t*sustr/500 and v = zeta = 0.
remora.prob_name = DoubleGyre
remora.solve3d = false
remora.max_step = 4
remora.prob_lo = 0 0 -500
remora.prob_hi = 16000 16000 0
remora.n_cell = 16 16 1
remora.is_periodic = 1 1 0
remora.fixed_dt = 10
remora.rdrag = 0
remora.use_coriolis = false
remora.horizontal_mixing_type = constant
remora.visc2 = 0
remora.rho0 = 1025
remora.g = 9.80665
remora.plot_file = plt
remora.plot_int = 4
remora.plot_vars_3d = temp salt x_velocity y_velocity z_velocity
remora.plot_vars_2d = zeta ubar vbar
remora.check_file = chk
remora.check_int = 2
amrex.init_snan = 1
amrex.fpe_trap_invalid = 1
