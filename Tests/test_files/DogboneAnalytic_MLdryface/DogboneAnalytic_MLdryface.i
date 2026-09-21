# ------------------  INPUTS TO MAIN PROGRAM  -------------------
# Asserts that a criterion on a face velocity is guarded by the face's mask, not the cell's.
#
# x_velocity is stored at a cell index but lives on that cell's low-x face, and vert_mean_3d
# multiplies it by msku(i,j) = mskr(i-1,j) * mskr(i,j). So a *water* cell whose i-1 neighbor is
# land holds an exact zero that is a mask artifact, not slack water. Guarding it by mskr(i,j)
# alone -- it is wet, so test it -- lets the jump between that zero and the real flow one cell
# over refine the coast, which is what the removed derefine criteria used to paper over.
#
# Nothing else in the suite sees this: MLcoastskip pins the guard for a cell-centered tracer,
# where the cell's own mask is the whole story, and deleting the per-field stencil leaves it,
# and every other case, green.
#
# The geometry puts a coast where there is real flow to contrast against. The free-surface bump
# drives x < 1100, so mask_x_lo/mask_x_hi move the land band inside it; at the shipped
# 2800/5600 the water by those coasts is nearly stagnant and the artifact zero is
# indistinguishable from its neighbors. regrid_int is positive because the masking that creates
# the artifact happens while stepping, not at init.
#
# The assertion is a cell count because that is what changes: the cell-only guard reaches the
# coastal columns and the face stencil does not, but both regions share a bounding box, so
# check_level_extent.sh cannot separate them. A cell count is summed over boxes and so survives
# load-balance chopping. This case is driven by the bump, not stationary, so there is no
# plt00000 comparison to make.
remora.prob_name = DogboneAnalytic

remora.max_step = 10

amrex.fpe_trap_invalid=1

# PROBLEM SIZE & GEOMETRY
remora.prob_lo     =   0.       0.  -10.
remora.prob_hi     =   8400.  750.       0.

remora.n_cell           =  42 15 16

# Pinned because this case asserts grid structure: each of these moves the assertion without
# any tagging having changed. main.cpp sets n_error_buf and blocking_factor imperatively rather
# than taking the AMReX defaults, and Inputs.rst documents a different n_error_buf than main.cpp
# sets, so a commit reconciling the two would otherwise turn these cases red for an unrelated
# reason. grid_eff matches the AMReX default, so pinning it is a no-op today.
amr.n_error_buf     = 0
amr.blocking_factor = 1
amr.grid_eff        = 0.7

amr.blocking_factor_z = 16

amr.max_grid_size_z = 1024

remora.is_periodic = 0 0 0

remora.bc.xlo.type = "slipwall"
remora.bc.xhi.type = "slipwall"
remora.bc.ylo.type = "slipwall"
remora.bc.yhi.type = "slipwall"

# TIME STEP CONTROL
remora.fixed_dt            = 2.0 # Timestep size (seconds)
remora.ndtfast = 20

# REFINEMENT / REGRIDDING
amr.max_level       = 1       # maximum level number allowed
amr.ref_ratio_vect = 3 3  1
# The artifact zero is created by the masking applied while stepping, so tagging has to happen
# after a step and not only inside InitFromScratch.
amr.regrid_int      = 2

# DIAGNOSTICS & VERBOSITY
remora.sum_interval  = 1
remora.v             = 0

# CHECKPOINT FILES
remora.check_file      = chk

# PLOTFILES
remora.plot_file     = plt
remora.plot_int      = 10
remora.plot_vars_3d  = salt temp x_velocity y_velocity z_velocity
remora.plot_vars_2d  = mask_rho   # a rho2d sidecar, for eyeballing a failure; fcompare ignores it
remora.plotfile_type = amrex

# SOLVER CHOICE
remora.tracer_horizontal_advection_scheme = "upstream3"

remora.Akt_bak = 1e-6
remora.Akv_bak = 1e-5

remora.use_coriolis  = false

remora.theta_s = 0.0
remora.theta_b = 0.0
remora.tcline = 1e16

remora.bottom_stress_type = "quadratic"
remora.rdrag2 = 3.0e-3

remora.mask_type = "analytic"

remora.init_l1ad_h =false
remora.init_l1ad_T =false

remora.init_l0int_h =true
remora.init_l0int_T =true

remora.init_ana_h = false
remora.init_ana_T = false

# PROBLEM PARAMETERS
remora.R0    = 1027.0
remora.S0    = 35.0
remora.T0    = 10.0
remora.Tcoef = 1.7e-4
remora.Scoef = 7.6e-4
remora.rho0  = 1025.0

remora.ic_type       = "analytic"

# The bump is the point here: it is what puts flow next to the coast.
remora.prob.zeta_bump = true
remora.prob.temp_west = 10.0

# Land band moved west, into the region the bump drives
remora.prob.mask_x_lo = 200.0
remora.prob.mask_x_hi = 600.0
remora.prob.mask_y_lo = 250.0
remora.prob.mask_y_hi = 500.0

# Above the real adjacent differences in the flow, below the jump from the artifact zero to
# the flow beside it, so the dry faces are the only thing that could tag here.
#
# Only the ratio of the field's differences to this threshold matters, so the assertion below
# holds over a band of thresholds rather than at a point. Measured on this build: 6480 cells
# for anything in roughly [0.019, 0.032], stepping to 8640 around 0.018 and 10800 by 0.016
# below the band, and to 5760 by 0.033 above it. 0.025 sits near the middle, about 30% from
# either edge, so a physics change would have to move the near-coast velocity differences by
# that much before this case turned red for a reason unrelated to tagging.
#
# If it does turn red, check that before suspecting the guard: this case has no gold file, so
# the cell count is its only signal and a changed flow reports as a tagging failure. Comparing
# plt00010 against a known-good run will tell you which it is.
remora.refinement_indicators = dryface
remora.dryface.max_level = 1
remora.dryface.adjacent_difference_greater = 0.025
remora.dryface.field_name = x_velocity

remora.coupling_type = "TwoWay"
