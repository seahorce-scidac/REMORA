# A uniform current of (0.5, 0.25) m/s with bulk fluxes. The test runs one step twice:
# with the wind taken relative to the current, under a wind of (5.5, 1.25) m/s, and with
# the absolute wind, under (5, 1) m/s. On the first step both see the same current, so
# both see the same relative wind, and the two must agree.
remora.max_step = 1
remora.plot_int = 1

remora.prob.u_0 = 0.5
remora.prob.v_0 = 0.25

remora.bulk_fluxes = true
