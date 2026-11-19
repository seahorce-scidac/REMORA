# Checkpoint step 1 precedes the layout change at time 2, before step 3.
remora.n_cell = 20 24 12
remora.fixed_dt = 1
remora.max_step = 3
remora.plot_int = 3
remora.check_int = 1
remora.do_substep = 0
amr.regrid_int = 1
amr.blocking_factor = 1
amr.n_error_buf = 0
remora.refinement_indicators = seam expanded
remora.expanded.max_level = 1
remora.expanded.in_box_lo = 0. 50000. -50.
remora.expanded.in_box_hi = 100000. 250000. 0.
remora.expanded.start_time = 2
remora.expanded.end_time = 3
