# As Advection_wind_minus_current, with half the current subtracted: the relative wind
# under (5.25, 1.125) m/s must match the absolute wind under (5, 1) m/s. The factor is
# set on the command line of the relative-wind run only.
remora.max_step = 1
remora.plot_int = 1

remora.prob.u_0 = 0.5
remora.prob.v_0 = 0.25

remora.bulk_fluxes = true
