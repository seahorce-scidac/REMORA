# Uniform inertial rotation changes the CFL timestep on every step.
remora.solve3d = false
remora.prob_hi = 16000 16000 0
remora.n_cell = 16 16 1
remora.use_coriolis = true
remora.coriolis_f0 = 0.01
remora.coriolis_beta = 0
remora.prob.v_0 = 0
remora.horizontal_mixing_type = constant
remora.visc2 = 0
remora.rdrag = 0
remora.fixed_dt = -1
remora.cfl = 0.5
remora.max_step = 4
remora.plot_int = 4
amrex.init_snan = 1
