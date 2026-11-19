# An embedding application supplies native AMReX geometry.
amr.n_cell = 8 8 8
geometry.prob_lo = 0 0 -10
geometry.prob_hi = 800 800 0
geometry.is_periodic = 1 1 0
amr.refine_grid_layout_z = 0
amr.max_grid_size_z = 8
amr.n_proper = 2
amr.max_grid_size = 4
