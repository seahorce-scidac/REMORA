# Compare the sheared profile with its discrete thickness-weighted depth mean.
# The initial diagnostic velocity must contain this mean at every depth.
remora.solve3d = false
remora.max_step = 0
remora.n_cell = 16 32 4
remora.plot_int = 1
remora.plot_vars_2d = zeta ubar vbar
remora.check_int = -1
amrex.init_snan = 1
remora.theta_s = 3
remora.theta_b = 0
remora.tcline = 25
