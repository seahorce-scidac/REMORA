GulfRefinementTest
==================

Developer lanes for ``remora.hires_grid_level`` and ``remora.hires_init_level``: bathymetry and
grid metrics (grid) or the initial state (init) specified on a refined level and averaged down to
level 0. These are the only lanes that exercise the NetCDF side of that feature -- the CI lanes
in ``Tests/test_files/`` cover the analytic side, and high-resolution *initialization* is
NetCDF-only (analytic aborts by design).

Requires a NetCDF build (``-DREMORA_ENABLE_PNETCDF=ON``, or ``USE_PNETCDF=TRUE`` for GNUmake).

The NetCDF files are not tracked in git (``.gitignore`` covers ``Exec/*/*.nc``). They come from
https://github.com/seahorce-scidac/REMORA-data.

Lanes
-----

``inputs_outflow``
    Gulf Stream toy grid, ``prob_name = IdealMiniGrid``, outflow on all four sides, no boundary
    file. ``n_cell = 31 9 100``, ``ref_ratio 3``, with **both** hires levels at 1: grid data from
    ``GS_toyGrd_lev1_v1.1.nc`` and initial data from ``GS_toyIni_lev1_v1.3.nc``. This is the lane
    that exercises the ``pm``/``pn`` average-down and rescaling, which the analytic lanes cannot
    reach (with analytic initialization ``set_grid_scale`` derives the metrics from the geometry
    instead).

``inputs_cf_orlanski``
    Gulf 5 km grid with Chapman/Flather/Orlanski boundary conditions from
    ``5km_bc_gulf_test_classic64.nc``. **Hires is off in this lane** -- see "Known data problems".

``inputs_synth_hires``
    Synthetic files from ``Tests/tools/make_hires_test_data.py`` (~200 KB, 64-bit offset classic),
    whose metrics ARE self-consistent -- unlike the GS toy level-1 grid. This is the lane to trust
    for numbers. Generate the data, then run::

        python3 ../../Tests/tools/make_hires_test_data.py --dir . --n-cell 10 10 4 --ref-ratio 3
        mpiexec -n 2 ./remora_exec inputs_synth_hires remora.max_step=0

    and check ``VOLUME`` at ``TIME = 0`` against the independent reference::

        python3 ../../Tests/reference/hires_volume_reference.py \
                --grid hires_grd.nc --init hires_ini.nc --ref-ratio 3

    With the default generator settings both give ``9522606688.29065``; measured agreement is
    4e-16 relative. That pins the average-down arithmetic: ``h`` and ``zeta`` averaged as plain
    r x r means, ``pm``/``pn`` averaged and then divided by ``r``, and the fine-to-coarse block
    alignment with the grow rings excluded. Because ``sum_k Hz == zeta + h`` exactly at t = 0,
    ``VOLUME`` is a direct readout of all four.

Grow-cell rule
--------------

A high-resolution file has to carry the refined domain plus ``cum_ref_ratios[hires_*_level]``
grow rings on every side::

    xi_rho  >= n_cell_x * r + 2r
    eta_rho >= n_cell_y * r + 2r

The GS toy set satisfies this exactly at ``r = 3``: level 0 is ``xi_rho = 33, eta_rho = 11``
(interior ``31 x 9``, matching ``n_cell``), level 1 is ``99 x 33``, and ``31*3 + 6 = 99``,
``9*3 + 6 = 33``. REMORA now checks this up front and aborts naming the file, the found size and
the required size; it used to read an undersized file cleanly and leave the outer rings unwritten.

Known data problems
-------------------

**The Gulf 1 km file is too small.** ``inputs_cf_orlanski`` has ``n_cell = 19 19`` with
``ref_ratio 5``, which needs ``19*5 + 2*5 = 105`` points per side, and
``1km_gulf_test_classic64.nc`` is ``101 x 101`` -- two rings short. The lane therefore keeps
``remora.hires_grid_level = -1`` until the file is regenerated at ``105 x 105``. (This defect was
latent: the input previously spelled the parameter ``remora.nc_hires_grid_level``, which ParmParse
silently ignores, so the lane looked like it exercised the feature and never opened the file.)

**``GS_toyGrd_lev1_v1.1.nc`` carries unrefined metrics.** It covers the same lon/lat extent as
``GS_toyGrd_lev0_v1.3.nc`` with 3x the cells per side, so its ``pm``/``pn`` should be 3x the
level-0 values -- but they are identical to them. Since the average-down divides ``pm`` by the
refinement ratio (correctly: coarsening multiplies ``dx`` by ``r``), the level-0 metrics that come
out of the hires path are 3x too small and cell areas 9x too large. Measured at ``max_step = 0``
with ``remora.sum_precision = 15``::

    hires_grid_level=1 hires_init_level=1   VOLUME = 2.72620210891816e+15
    hires off (level-0 files only)          VOLUME = 3.02911427291809e+14
    ratio                                   9.0000  == r^2

and directly from the files, summing ``h/(pm*pn)`` over each interior::

    GS_toyGrd_lev0_v1.3.nc   area 7.87397e+10   volume 3.02921e+14
    GS_toyGrd_lev1_v1.1.nc   area 7.08792e+11   volume 2.72618e+15   <- 9x, same pm as lev0

So the ``r^2`` discrepancy is in the input file, not in the average-down: the code reads ``pm``
from the file, averages it down and rescales by ``r`` as it should. Note this file is the only one
in the set tagged ``v1.1`` while its companions are ``v1.3``. **Regenerate it with refined metrics
before treating this lane's numbers as physically meaningful**; until then it is a plumbing test
(does the read, average-down and rescale run and stay finite), not a validation.

Running
-------

::

    mpiexec -n 2 ./remora_exec inputs_outflow remora.max_step=2 \
            remora.v=1 remora.sum_interval=1 remora.sum_precision=15

To confirm the hires path is actually doing something, run the same input with
``remora.hires_grid_level=-1 remora.hires_init_level=-1`` and compare ``VOLUME`` at ``TIME`` 0.
Identical values would mean the hires files were never read.

Stale files in this directory
-----------------------------

``tmp_build_dir/`` (323 MB), ``CMakeLists.txt``, ``GNUmakefile``, ``Make.package`` and
``AMReX_buildInfo.cpp`` were copied in from a pre-unification tree and are dead: they reference a
``prob.H``/``prob.cpp`` that no longer exists anywhere in the repo, ``CMakeLists.txt`` sets
``remora_exe_name`` to ``idealminigrid``, and nothing under ``Exec/`` is added as a CMake
subdirectory any more -- the build produces a single ``remora_exec`` and the problem is selected at
run time by ``remora.prob_name``. They can be removed::

    rm -rf Exec/GulfRefinementTest/tmp_build_dir
    rm -f  Exec/GulfRefinementTest/{CMakeLists.txt,GNUmakefile,Make.package,AMReX_buildInfo.cpp}
