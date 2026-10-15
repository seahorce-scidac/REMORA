
 .. _RegressionTests:

Regression Tests
================

There are currently 63 regression tests which are run as part of every PR on a range of architectures:

* 22 that compare a plotfile against a stored reference;
* 6 that compare two runs of the model against each other, or a run against its own initial state, and
  need no reference;
* 5 that assert values against an external reference;
* 5 that assert the structure of the refined grid;
* 8 that assert how far an integrated quantity drifts over a run;
* 8 that assert a run reaches a behavior it reports in its log;
* 9 that assert a misconfigured input fails for the stated reason.

There is also one unit test, ``DateClock``. The CI tests use cmake and are based on the version
of AMReX in the REMORA submodule. This suite can be run following the
instructions in :ref:`Testing<Testing>`, which also describes the different kinds of assertion and what
each is for. All of them run under MPI on two ranks, and are checked to be rank-invariant.

In addition there is a suite of more extensive nightly tests that use GNUMake and use the current
development branch of AMReX.

Results from the nightly CPU tests can be found here: `CPU tests`_

Results from the nightly GPU tests can be found here: `GPU tests`_

.. _`CPU tests`: https://ccse.lbl.gov/pub/RegressionTesting1/REMORA

.. _`GPU tests`: https://ccse.lbl.gov/pub/GpuRegressionTesting/REMORA

Continuous Integration (CI) Tests
---------------------------------

The following problems are currently tested in the CI. More details about the problems underlying these tests are given in :ref:`sec:verification`.

+-------------------------+----------+-----------+----------+---------------------------------+
| Test                    | nx ny nz | xbc       | ybc      | Other                           |
+=========================+==========+===========+==========+=================================+
| Advection               | 81 81 16 | Periodic  | Periodic |                                 |
+-------------------------+----------+-----------+----------+---------------------------------+
| Advection_ML            | 80 80 16 | Periodic  | Periodic | multilevel                      |
+-------------------------+----------+-----------+----------+---------------------------------+
| BoundaryLayer           | 39 4 30  | Radiation | Periodic | Coriolis                        |
|                         |          |           |          |                                 |
|                         |          | / outflow |          | GLS mixing scheme               |
|                         |          |           |          |                                 |
|                         |          |           |          | non-flat bathymetry             |
+-------------------------+----------+-----------+----------+---------------------------------+
| Channel_Test            | 20 60 50 | Periodic  | SlipWall | Coriolis                        |
|                         |          |           |          |                                 |
|                         |          |           |          | GLS mixing scheme               |
|                         |          |           |          |                                 |
|                         |          |           |          | non-flat bathymetry             |
|                         |          |           |          |                                 |
|                         |          |           |          | quadratic bottom stress         |
|                         |          |           |          |                                 |
|                         |          |           |          | bulk fluxes                     |
|                         |          |           |          |                                 |
|                         |          |           |          | cloud cover                     |
|                         |          |           |          |                                 |
|                         |          |           |          | evaporation/precipitation with  |
|                         |          |           |          |                                 |
|                         |          |           |          | sea surface height correction   |
+-------------------------+----------+-----------+----------+---------------------------------+
| Channel_Test_ML_seam    | 20 60 50 | Periodic  | SlipWall | Coriolis                        |
|                         |          |           |          |                                 |
|                         |          |           |          | static multilevel, refined      |
|                         |          |           |          |                                 |
|                         |          |           |          | patch on the periodic seam      |
+-------------------------+----------+-----------+----------+---------------------------------+
| DogboneAnalytic         | 42 15 16 | SlipWall  | SlipWall | quadratic bottom stress         |
|                         |          |           |          |                                 |
|                         |          |           |          | land-sea masking                |
+-------------------------+----------+-----------+----------+---------------------------------+
| DogboneAnalytic_MLquad  | 42 15 16 | SlipWall  | SlipWall | quadratic bottom stress         |
|                         |          |           |          |                                 |
|                         |          |           |          | land-sea masking                |
|                         |          |           |          |                                 |
|                         |          |           |          | static multilevel               |
+-------------------------+----------+-----------+----------+---------------------------------+
| DogboneAnalytic_MLvel   | 42 15 16 | SlipWall  | SlipWall | quadratic bottom stress         |
|                         |          |           |          |                                 |
|                         |          |           |          | land-sea masking                |
|                         |          |           |          |                                 |
|                         |          |           |          | dynamic multilevel              |
+-------------------------+----------+-----------+----------+---------------------------------+
| DoubleGyre              | 54 108 4 | SlipWall  | SlipWall | Coriolis                        |
+-------------------------+----------+-----------+----------+---------------------------------+
| DoublyPeriodic          | 41 80 16 | Periodic  | Periodic | Coriolis                        |
|                         |          |           |          |                                 |
|                         |          |           |          | non-flat bathymetry             |
+-------------------------+----------+-----------+----------+---------------------------------+
| Seamount                | 49 48 13 | Periodic  | Periodic | Coriolis                        |
+-------------------------+----------+-----------+----------+---------------------------------+
| Upwelling               | 41 80 16 | Periodic  | SlipWall | Coriolis                        |
|                         |          |           |          |                                 |
|                         |          |           |          | non-flat bathymetry             |
+-------------------------+----------+-----------+----------+---------------------------------+
| Upwelling_Fennel        | 41 80 16 | Periodic  | SlipWall | Coriolis                        |
|                         |          |           |          |                                 |
|                         |          |           |          | non-flat bathymetry             |
|                         |          |           |          |                                 |
|                         |          |           |          | bulk fluxes                     |
|                         |          |           |          |                                 |
|                         |          |           |          | Fennel biology, analytic IC     |
|                         |          |           |          |                                 |
|                         |          |           |          | carbon, oxygen, river DON       |
|                         |          |           |          |                                 |
|                         |          |           |          | non-conservative alkalinity     |
|                         |          |           |          |                                 |
|                         |          |           |          | dated atmospheric pCO2          |
|                         |          |           |          |                                 |
|                         |          |           |          | Wanninkhof (2014) gas transfer  |
+-------------------------+----------+-----------+----------+---------------------------------+
| Upwelling_GLS           | 41 80 16 | Periodic  | SlipWall | Coriolis                        |
|                         |          |           |          |                                 |
|                         |          |           |          | non-flat bathymetry             |
|                         |          |           |          |                                 |
|                         |          |           |          | GLS mixing scheme               |
+-------------------------+----------+-----------+----------+---------------------------------+
| Upwelling_NLEOS         | 41 80 16 | Periodic  | SlipWall | Coriolis                        |
|                         |          |           |          |                                 |
|                         |          |           |          | non-flat bathymetry             |
|                         |          |           |          |                                 |
|                         |          |           |          | nonlinear equation of state     |
+-------------------------+----------+-----------+----------+---------------------------------+
| Upwelling_logdrag       | 41 80 16 | Periodic  | SlipWall | Coriolis                        |
|                         |          |           |          |                                 |
|                         |          |           |          | non-flat bathymetry             |
|                         |          |           |          |                                 |
|                         |          |           |          | logarithmic bottom stress       |
+-------------------------+----------+-----------+----------+---------------------------------+
| Upwelling_qdrag         | 41 80 16 | Periodic  | SlipWall | Coriolis                        |
|                         |          |           |          |                                 |
|                         |          |           |          | non-flat bathymetry             |
|                         |          |           |          |                                 |
|                         |          |           |          | quadratic bottom stress         |
+-------------------------+----------+-----------+----------+---------------------------------+

``Channel_Test_ML_seam`` places a refined patch whose low x edge lies on the periodic seam. It is the only
lane whose fine level touches a periodic boundary, and so the only one covering the periodic handling in
the mask construction, the ghost fill of the fine level, and the ubar/vbar average-down. It is also the
only multi-level lane with Coriolis on. Its reference was produced on one rank with one box per level;
under MPI the level-1 box splits, so the two-rank run also covers the ghost refresh at the internal box
boundary after the coarse-fine barotropic boundary conditions are set.

.. _sec:ci-hires-tests:

High-resolution initialization tests
------------------------------------

These cover ``remora.hires_grid_level`` and ``remora.hires_init_level``: bathymetry or the initial state
evaluated on a refined level and averaged down to level 0 (see :ref:`Inputs<sec:Inputs>`). They are
analytic, so they need no NetCDF build and no external data. The NetCDF reader behind
``remora.hires_init_level`` is covered by the developer lanes in ``Exec/GulfRefinementTest`` instead,
described in that directory's ``README.rst``.

The lanes come in three flavors, because "the run completed" is not evidence that an average-down
happened, let alone that it was right.

*Transparency.* ``ChannelTest`` sets ``h = 50`` and ``DogboneAnalytic`` sets ``h = 10``, both with
``setVal``, covering every grow cell. The mean over any refined block is then the same constant, so the
hires path must reproduce the corresponding non-hires baseline exactly. These lanes compare against the
existing ``Channel_Test`` and ``DogboneAnalytic_MLvel`` reference files rather than storing new ones.
They deliberately *cannot* detect a hires flag that is being ignored -- that is the next flavor's job.

*Discrimination.* ``Seamount``'s bathymetry is a Gaussian in physical position, so the mean of the
refined samples in a coarse cell is not the coarse midpoint sample. ``Seamount_hires`` must therefore
both match its own reference and **differ** from the plain ``Seamount`` one; the difference is about
16 m at the summit, some 12 orders of magnitude above the comparison tolerance. Without that second
clause a feature that silently stops taking effect still matches its own snapshot.

*Reference values.* The three ``Seamount_hires_init`` lanes assert level-0 temperature extrema computed by
``Tests/reference/hires_init_analytic_reference.py``, which rebuilds the Seamount bathymetry, vertical
grid, and temperature profile in NumPy on the refined grid and averages them down. It reproduces the
whole level-0 temperature field of all three lanes to about 1e-14, against a gap of about 7e-2 to the
plain level-0 evaluation. The asserted maxima therefore fail for a run that ignores
``hires_init_level`` (by 1.4e-4 against a 1e-8 tolerance), and they are not snapshots of the code under
test. The three lanes take the bathymetry at the initialization level from a different source each:
evaluated there, evaluated there although grid data exists on a coarser level, and averaged down from
grid data on a finer level.

*Misconfiguration.* Four lanes assert that a bad combination fails with the message that names it, rather
than segfaulting, writing out of bounds, or running on data it quietly ignored.

*Partial coverage.* The two ``DogboneAnalytic_MLmask`` tests are the only ones in which a coarse cell is
part land and part water, which is the case the mask-weighted average-down exists for. Producing such a
cell in the Dogbone case takes two things: ``hires_grid_level``, which resolves the coastline on the
refined level instead of injecting it from level 0, and a static box that places the fine grids over
the coast. The box is static because these cases are at rest, so a field-based indicator such as
``DogboneAnalytic_MLvel``'s ``x_velocity > 0.05`` would tag nothing anywhere; it also fixes which coarse
cells come out partly wet, and with ``amr.regrid_int = -1`` they cannot move between the two plotfiles.
Being at rest is also why they need no reference file of their own -- ``plt00010`` must equal ``plt00000``
to 1e-14. Under a plain arithmetic average-down the zeroed fine land cells drag the partly covered coarse
cells off that initial value and stationarity breaks by about 1e-2 in salinity and velocity, so a
lost mask weighting fails the comparison rather than quietly agreeing with its own snapshot.

+---------------------------------+----------+---------------------------------------------------------+
| Test                            | nx ny nz | What it asserts                                         |
+=================================+==========+=========================================================+
| Channel_Test_hires              | 20 60 50 | constant bathymetry averages down exactly (ratio 3),    |
|                                 |          |                                                         |
|                                 |          | compared against the Channel_Test reference             |
+---------------------------------+----------+---------------------------------------------------------+
| DogboneAnalytic_MLhires         | 42 15 16 | the same, with a live refined level and regridding,     |
|                                 |          |                                                         |
|                                 |          | so it also covers the regrid paths and masking          |
+---------------------------------+----------+---------------------------------------------------------+
| Seamount_hires                  | 49 48 13 | a Gaussian bathymetry averages down to something        |
|                                 |          |                                                         |
|                                 |          | different from the direct level-0 evaluation (ratio 3)  |
+---------------------------------+----------+---------------------------------------------------------+
| Seamount_hires_r4               | 49 48 13 | the same at ratio 4, wider than the standard grow       |
|                                 |          |                                                         |
|                                 |          | ring, and differing from the ratio-3 result             |
+---------------------------------+----------+---------------------------------------------------------+
| Upwelling_Fennel_hires_init     | 41 80 16 | biology coexists with an averaged-down bathymetry:      |
|                                 |          |                                                         |
|                                 |          | the six constant Fennel tracers and the zero dye are    |
|                                 |          |                                                         |
|                                 |          | asserted by value, catching a component-offset error    |
+---------------------------------+----------+---------------------------------------------------------+
| Seamount_hires_init             | 49 48 13 | the initial state evaluated on level 1 and averaged     |
|                                 |          |                                                         |
|                                 |          | down matches the reference extrema (ratio 3)            |
+---------------------------------+----------+---------------------------------------------------------+
| Seamount_hires_init_above_grid  | 49 48 13 | the same with initialization on level 2 above grid data |
|                                 |          |                                                         |
|                                 |          | on level 1, over the analytic level-2 bathymetry (6)    |
+---------------------------------+----------+---------------------------------------------------------+
| Seamount_hires_init_below_grid  | 49 48 13 | the same with initialization on level 1 below grid data |
|                                 |          |                                                         |
|                                 |          | on level 2, over the averaged-down bathymetry           |
+---------------------------------+----------+---------------------------------------------------------+
| DogboneAnalytic_hires_init\_    | 42 15 16 | zeta averaged down from level 1: VOLUME at t = 0 equals |
| zeta                            |          |                                                         |
|                                 |          | its closed-form value, 63566666.667                     |
+---------------------------------+----------+---------------------------------------------------------+
| Upwelling_Fennel_hires_init_ic  | 41 80 16 | biology evaluated on level 1 and averaged down: the     |
|                                 |          |                                                         |
|                                 |          | Fennel constants survive exactly, and the max           |
|                                 |          |                                                         |
|                                 |          | temperature separates it from the lane above            |
+---------------------------------+----------+---------------------------------------------------------+
| Seamount_hires_init_gridscale\_ | 49 48 13 | analytic high-resolution initialization with a          |
| abort                           |          |                                                         |
|                                 |          | non-constant grid scale is rejected                     |
+---------------------------------+----------+---------------------------------------------------------+
| Seamount_hires_grid_max_abort   | 49 48 13 | a hires level above ``amr.max_level`` is rejected       |
+---------------------------------+----------+---------------------------------------------------------+
| Seamount_hires_init_max_abort   | 49 48 13 | the same for the initialization level                   |
+---------------------------------+----------+---------------------------------------------------------+
| Seamount_hires_grid_zero_abort  | 49 48 13 | a hires level of 0 is rejected rather than              |
|                                 |          |                                                         |
|                                 |          | dereferencing an array only allocated above level 0     |
+---------------------------------+----------+---------------------------------------------------------+
| DogboneAnalytic_MLmask          | 42 15 16 | coarse cells that are 6 water of 9 keep their exact     |
|                                 |          |                                                         |
|                                 |          | values, so a run at rest stays at rest (ratio 3)        |
+---------------------------------+----------+---------------------------------------------------------+
| DogboneAnalytic_MLmask_rr2      | 42 15 16 | the same at ratio 2, where the blocks are half water    |
+---------------------------------+----------+---------------------------------------------------------+
| DogboneAnalytic_MLmask_hires\_  | 42 15 16 | the same with the initial state also averaged down      |
| init                            |          |                                                         |
|                                 |          | from level 1 through the partly wet cells               |
+---------------------------------+----------+---------------------------------------------------------+

``Upwelling_Fennel_hires_init`` runs to zero steps: it writes the initial plotfile and exits, which takes
about a second and keeps the time stepper out of an assertion about initialization. The six Fennel
constants come from the analytic profile in ``Source/Prob/REMORA_InitAnalyticBiology_BioToy.H``, and a
constant is unchanged by averaging, so the expected level-0 values are known outright rather than
snapshotted. Asserting the dye is zero alongside them is what catches a regression in the
``Bio_comp = Tracer_comp + remora.nscalar`` offset, which would shift a biology tracer into the dye slot.

The two ``DogboneAnalytic_MLmask`` lanes also assert where level 1 is. The static box spans the full
width of the domain in y, 0 to 750, and the coasts run through it, so the whole box must be refined,
land included, rather than clipped back to the water (which once gave y in [250, 500]). The ratio-2 lane
covers the same fix on the other refinement-ratio path: there the grids go from three clipped boxes to
two full-box grids, which a stationarity comparison alone cannot see.

.. _sec:ci-tagging-tests:

Refinement tagging tests
------------------------

Grid structure is the only place a tagging bug shows: a criterion that has stopped tagging, or one that
has started tagging the coastline, still writes a self-consistent plotfile that a snapshot comparison
accepts. These lanes assert the grid itself instead -- which levels exist, or how many cells a level
covers -- using ``adjacent_difference_greater`` criteria on the Dogbone case at ratio 3.

+---------------------------------+----------+---------------------------------------------------------+
| Test                            | nx ny nz | What it asserts                                         |
+=================================+==========+=========================================================+
| DogboneAnalytic_MLcoastskip     | 42 15 16 | a gradient criterion on ``temp`` does not see the jump  |
|                                 |          |                                                         |
|                                 |          | to the zeroed land cells: no level 1 is built, and the  |
|                                 |          |                                                         |
|                                 |          | run stays at rest                                       |
+---------------------------------+----------+---------------------------------------------------------+
| DogboneAnalytic_MLcoasttag      | 42 15 16 | the same criterion on ``mask`` does see the coast and   |
|                                 |          |                                                         |
|                                 |          | builds level 1                                          |
+---------------------------------+----------+---------------------------------------------------------+
| DogboneAnalytic_MLdryface       | 42 15 16 | a criterion on ``x_velocity`` is guarded by the masks   |
|                                 |          |                                                         |
|                                 |          | of both cells sharing the face, not by the cell's own:  |
|                                 |          |                                                         |
|                                 |          | level 1 covers exactly 6480 cells                       |
+---------------------------------+----------+---------------------------------------------------------+

``MLcoastskip`` and ``MLcoasttag`` are a pair: ``MLcoastskip`` catches a guard that has stopped skipping
the coast, and ``MLcoasttag`` catches a lost ``mask`` exemption. Neither catches a guard that tags nothing
at all. The exemption is a null mask, so ``MLcoasttag`` runs AMReX's unguarded test, and ``MLcoastskip``
expects no level 1 anyway. That case is caught by ``DogboneAnalytic_MLvel`` and ``_MLhires``, whose
``x_velocity`` criteria go through the guard; ``Advection_ML`` does not, since an unmasked run skips it.

``MLcoastskip`` and ``MLdryface`` set ``amr.regrid_int = 2``. The land values they contrast against are
created by masking during stepping -- at initialization ``temp`` is 10 on land as well as water -- so
tagging only at initialization would let both pass with the guard removed.

``MLdryface`` covers a face-staggered field, where the cell's own mask is not the whole story. A wet cell
whose low-x neighbor is land holds an exact zero in ``x_velocity`` that is a mask artifact, not slack
water. Guarding it by that cell's mask alone lets the jump between the zero and the real flow next to it
refine the coast, which refines 12960 cells instead of 6480. Both regions share a bounding box, so the
assertion is a cell count, summed over boxes so it does not depend on load balancing. The land band is
moved into the region the free-surface bump drives, so there is real flow beside the coast to contrast
against, and this case evolves, so there is no stationarity comparison. Its threshold, 0.025, sits near
the middle of the band that gives 6480 cells (roughly [0.019, 0.032]; 8640 or more below it, 5760 above),
so the flow has to change by about 30% before the count moves. With no gold file the count is the only
signal, so a failure here may be a change in the flow rather than in tagging.

The ``vorticity`` stencil, a 3x3 block, has no test: the only masked problem is near-irrotational, so any
assertion on it would be pinned to roundoff.

.. _sec:ci-subcycle-tests:

Time subcycling tests
---------------------

These cover ``remora.do_substep``, under which each refined level advances with its own time step. The
answers on the default path are pinned by the ``Advection_ML`` and ``DogboneAnalytic_ML`` reference lanes;
these lanes assert what those cannot -- that subcycling actually engages, that the subcycling and lockstep
drivers agree where they should, and that the assumptions underneath the coarse-fine treatment hold.
Unless noted otherwise, level 1 runs at half the level-0 step (``fixed_dt = 100``, ratio 2).

+---------------------------------+----------+---------------------------------------------------------+
| Test                            | nx ny nz | What it asserts                                         |
+=================================+==========+=========================================================+
| Advection_ML_subcycle           | 80 80 16 | level 1 is advanced with ``dt = 50``                    |
+---------------------------------+----------+---------------------------------------------------------+
| Advection_ML_subcycle\_         | 80 80 16 | at a time-step ratio of 1, the subcycling driver        |
| identity                        |          |                                                         |
|                                 |          | reproduces the lockstep driver                          |
+---------------------------------+----------+---------------------------------------------------------+
| DogboneAnalytic_ML_subcycle     | 42 15 16 | subcycling engages at ratio 3 over varying bathymetry,  |
|                                 |          |                                                         |
|                                 |          | three level-1 steps of 0.6667 s per level-0 step        |
+---------------------------------+----------+---------------------------------------------------------+
| Advection_ML_cf_metrics         | 80 80 16 | the fine cell edges along the coarse-fine interface sum |
|                                 |          |                                                         |
|                                 |          | to the coarse edge, to 1e-14 (ratio 2)                  |
+---------------------------------+----------+---------------------------------------------------------+
| DogboneAnalytic_ML_cf_metrics   | 42 15 16 | the same at ratio 3, the grid the volume lanes below    |
|                                 |          |                                                         |
|                                 |          | are measured on                                         |
+---------------------------------+----------+---------------------------------------------------------+
| BoundaryLayer_ML_cf_metrics     | 39 4 30  | the same over non-uniform ``pm`` and ``pn``, where it   |
|                                 |          |                                                         |
|                                 |          | does not hold trivially                                 |
+---------------------------------+----------+---------------------------------------------------------+
| BoundaryLayer_ML_plot_z         | 39 4 30  | the node positions a viewer rebuilds from the plotfile  |
|                                 |          |                                                         |
|                                 |          | reproduce ``z_phys_nd`` on every level                  |
+---------------------------------+----------+---------------------------------------------------------+
| Advection_ML_do_substep_alias   | 80 80 16 | the old ``amr.do_substep`` spelling is still read, and  |
|                                 |          |                                                         |
|                                 |          | warns that it is deprecated                             |
+---------------------------------+----------+---------------------------------------------------------+
| Advection_ML_do_substep_both\_  | 80 80 16 | setting both ``remora.do_substep`` and                  |
| abort                           |          |                                                         |
|                                 |          | ``amr.do_substep`` is rejected                          |
+---------------------------------+----------+---------------------------------------------------------+

``Advection_ML_subcycle_identity`` is the load-bearing one: it separates a broken driver from the answer
changes subcycling legitimately makes. ``remora.do_reflux`` is off in its subcycled run, because
refluxing is a correction the lockstep driver does not apply at all, and leaving it on would compare a
feature rather than the drivers. ``Advection`` has a flat bottom, so the water depth matches across the
interface and the mass-flux boundary condition reduces to the interpolation it replaces; the Dogbone
lane, with varying bathymetry and ratio 3, is the one that exercises it.

The three ``cf_metrics`` lanes check an identity the conservative coarse-fine flux depends on: unless the
fine edges sum to the coarse edge, the mass flux imposed at the interface cannot be conservative. The
check runs under ``remora.check_cf_metrics`` and aborts past ``remora.check_cf_tol``, so reaching the
printed summary is the assertion. It measures 0 at ratio 2 and 1.4e-16 at ratio 3; every other
multi-level case has uniform ``pm`` and ``pn``, which is why the BoundaryLayer lane is needed. It is not
guaranteed on the NetCDF grid path, where a finer level's metrics are interpolated from the parent's.

``BoundaryLayer_ML_plot_z`` rebuilds each node's height the way a viewer does, from the plotfile header's
``prob_lo`` and ``dx`` plus the ``Nu_nd`` displacement, and aborts under ``remora.check_plot_z`` unless the
top and bottom node layers of every level reproduce ``z_phys_nd``. With ``expand_plotvars_to_unif_rr`` the
refined level's header ``dz`` is the native one divided by the ratio, and the refined patch's surface once
came out 32 m low because the nodal data was not refined in z to match.

.. _sec:ci-conservation-tests:

Conservation tests
------------------

``Advection`` is doubly periodic and ``DogboneAnalytic`` is closed by slip walls, so nothing can leave
either domain and the totals of tracer and volume have to hold. Each lane runs 20 steps and asserts the
relative drift of a total in the integrated-quantity log against a bound, either below it or, for a
control, above it. The bounds come from measurement.

+---------------------------------------+------------+------------------------------------------------+
| Test                                  | Quantity   | What it asserts                                |
+=======================================+============+================================================+
| Advection_ML_conservation             | tracer     | with refluxing, the drift is below 1e-10       |
+---------------------------------------+------------+------------------------------------------------+
| Advection_ML_conservation_control     | tracer     | without refluxing, the same run drifts above   |
|                                       |            |                                                |
|                                       |            | 1e-6                                           |
+---------------------------------------+------------+------------------------------------------------+
| Advection_3L_conservation             | tracer     | the same as the first, with three levels       |
+---------------------------------------+------------+------------------------------------------------+
| Advection_3L_conservation_control     | tracer     | the same as the second, with three levels      |
+---------------------------------------+------------+------------------------------------------------+
| Advection_conservation_baseline       | tracer,    | on one level, both drift below 1e-12           |
|                                       |            |                                                |
|                                       | volume     |                                                |
+---------------------------------------+------------+------------------------------------------------+
| DogboneAnalytic_ML_conservation       | volume     | the default subcycling path drifts below 1e-7  |
+---------------------------------------+------------+------------------------------------------------+
| DogboneAnalytic_ML_conservation_cons  | volume     | with the two ghost-band free-surface options   |
|                                       |            |                                                |
|                                       |            | off, the drift is below 1e-8                   |
+---------------------------------------+------------+------------------------------------------------+
| DogboneAnalytic_ML_conservation\_     | volume     | the lockstep driver drifts above 1e-7          |
| lockstep                              |            |                                                |
+---------------------------------------+------------+------------------------------------------------+

The controls are what make the bounds mean something: if one ever passes by conserving, its case has
stopped exercising the correction its partner is asserting, and the partner is proving nothing.
Refluxing takes the two-level tracer drift from 7.8e-5 to below what 12 digits resolve. The three-level
lanes exist because the correction at the level-1/level-2 interface has to be accumulated over several
level-1 steps before it is applied; a two-level case refluxes once per step of its only coarse level,
and so passes whether or not that accumulation is right. The single-level baseline is exact, so the AMR
bounds are measured against roundoff rather than against an unknown scheme error.

Volume is the integral of the free surface, so what costs it is how that surface is written in the fine
ghost band. By default REMORA writes it as ROMS's ``put_refine2d`` does -- every leapfrog record plus
``Zt_avg1`` (``remora.cf_fill_all_kcomp``), interpolated onto the child's sub-time
(``remora.cf_time_interp_zeta``) -- and neither preserves volume: the default path measures 2.9e-8, and
with both off the drift falls to 2.0e-9. The lockstep driver interpolates the parent's ubar at the
interface instead of imposing its time-averaged flux and drifts 2.0e-6, which is what separates the two
drivers in this case.

.. _sec:ci-ndtfast-tests:

Barotropic substep count tests
------------------------------

``remora.ndtfast`` is the number of barotropic steps per baroclinic step. It has no usable default, and a
zero divides by zero and sizes the barotropic weights to one element, which the 2D step then reads past
the end of. These lanes, all on ``Channel_Test``, cover it and the input names around it.

+---------------------------------+----------+---------------------------------------------------------+
| Test                            | nx ny nz | What it asserts                                         |
+=================================+==========+=========================================================+
| Channel_Test_ndtfast_alias      | 20 60 50 | ``remora.fixed_ndtfast_ratio`` gives the same answer as |
|                                 |          |                                                         |
|                                 |          | ``remora.ndtfast`` with the same count                  |
+---------------------------------+----------+---------------------------------------------------------+
| Channel_Test_ndtfast_unset\_    | 20 60 50 | a run that sets neither is rejected                     |
| abort                           |          |                                                         |
+---------------------------------+----------+---------------------------------------------------------+
| Channel_Test_fixed_fast_dt\_    | 20 60 50 | the removed ``remora.fixed_fast_dt`` is rejected rather |
| abort                           |          |                                                         |
|                                 |          | than silently ignored                                   |
+---------------------------------+----------+---------------------------------------------------------+
| Channel_Test_ndtfast_both\_     | 20 60 50 | setting both spellings is rejected                      |
| abort                           |          |                                                         |
+---------------------------------+----------+---------------------------------------------------------+

The alias comparison cannot pass for the wrong reason: a run that ignored the alias would be left without
a count and abort rather than agree.

.. _sec:ci-maskring-tests:

Wall masking tests
------------------

These place land on the boundary ring of rho points beside interior water, as in a ROMS closed basin
(``remora.prob.mask_ring``), on the ``DogboneAnalytic`` case with constant horizontal mixing so that the
ring is read.

+---------------------------------+----------+---------------------------------------------------------+
| Test                            | nx ny nz | What it asserts                                         |
+=================================+==========+=========================================================+
| DogboneAnalytic_maskring        | 42 15 16 | the wall boundary conditions put zero on that land, as  |
|                                 |          |                                                         |
|                                 |          | ROMS's ``t3dbc`` and ``zetabc`` do                      |
+---------------------------------+----------+---------------------------------------------------------+
| DogboneAnalytic_maskring\_      | 42 15 16 | writing a plotfile every step gives the same answer as  |
| plotint                         |          |                                                         |
|                                 |          | writing one at the end                                  |
+---------------------------------+----------+---------------------------------------------------------+
| DogboneAnalytic_maskring\_      | 42 15 16 | a run restarted from its step-5 checkpoint finishes     |
| restart                         |          |                                                         |
|                                 |          | where the continuous run does                           |
+---------------------------------+----------+---------------------------------------------------------+
| DogboneAnalytic_MLmaskring      | 42 15 16 | a refined patch against the walls, with its mask        |
|                                 |          |                                                         |
|                                 |          | injected from level 0, agrees with the same run with    |
|                                 |          |                                                         |
|                                 |          | ``hires_grid_level = 1``, so its grow cells beyond the  |
|                                 |          |                                                         |
|                                 |          | domain keep the ring                                    |
+---------------------------------+----------+---------------------------------------------------------+
| DogboneAnalytic_MLmaskring\_    | 42 15 16 | the same against a periodic x seam with land on its     |
| periodic                        |          |                                                         |
|                                 |          | far side, so those grow cells hold what lies across it  |
+---------------------------------+----------+---------------------------------------------------------+

Only the first needs a reference file.

.. _sec:ci-bc-tests:

Boundary condition tests
------------------------

These cover the two ways of specifying boundary conditions, per side and per variable (see
:ref:`Physical/Domain Boundary Conditions<sec:domainBCs>`). Both reach the same per-face state, so the
assertion is that they agree and neither needs a stored reference.

+---------------------------------+----------+---------------------------------------------------------+
| Test                            | nx ny nz | What it asserts                                         |
+=================================+==========+=========================================================+
| BC_per_variable                 | 54 108 4 | a per-variable ``West South East North`` list lands on  |
|                                 |          |                                                         |
|                                 |          | the same four faces the per-side keywords name, and     |
|                                 |          |                                                         |
|                                 |          | both routes read the same inflow values                 |
+---------------------------------+----------+---------------------------------------------------------+
| BC_inflow_no_value              | 54 108 4 | an inflow face with no inflow value for one of its      |
|                                 |          |                                                         |
|                                 |          | tracers is rejected, naming the input it wants          |
+---------------------------------+----------+---------------------------------------------------------+

``BC_per_variable`` runs both inputs in its test directory and compares the two plotfiles. All four sides
carry a different condition -- inflow, noslipwall, outflow, slipwall going West, South, East, North -- so
any other order changes the answer. A y-symmetric case cannot do this, so this lane is the one that targets a South/North
transposition.

The y pairing is discriminated by the *velocity* conditions specifically. For cell-centered tracers, and
for ``zeta`` and ``tke``, ``noslipwall`` and ``slipwall`` both translate to ``foextrap``, leaving South and
North indistinguishable in those fields; only ``u``/``v`` and ``ubar``/``vbar`` separate them, since
``no_slip_wall`` gives ``ext_dir`` on both components while ``slip_wall`` gives ``foextrap`` tangentially
and ``ext_dir`` normally. Substituting conditions that look equally distinct but coincide for the
velocities would blind this lane without changing its appearance.

The dye is a second readout: ``DoubleGyre`` initializes it to zero and the only dye enters through the
western inflow value, so a nonzero tracer field means that value was read, and an agreeing one means both
routes read it the same. Agreement alone cannot assert this -- two runs that both ignored the value would
agree at zero -- so the lane also carries a ``check_extrema.sh`` assertion on the dye's magnitude, with a
tolerance loose enough to be a claim about order of magnitude rather than a blessed digit string.

The inflow is driven by horizontal diffusion rather than advection. ``ubar`` and ``vbar`` take no inflow
value of their own, so the barotropic normal velocity at the face stays zero and the barotropic correction
holds the net flux near zero regardless of the 3D ``velocity`` entry: raising it tenfold moves the dye by
about 7%. Advective inflow is therefore not covered by either lane.

Nightly Regression Tests on CPU
-------------------------------
And the following are currently tested nighly on CPU.

Based on :ref:`Advection<advection>`:

+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Test                                   | nx ny nz     | xbc              | ybc               | Other                            |
+========================================+==============+==================+===================+==================================+
| Advection-1grid-xy                     | 81 81 16     | Periodic         | Periodic          |                                  |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Advection-1grid-xy-ML                  | 80 80 16     | Periodic         | Periodic          | multilevel                       |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Advection-OMP-xy                       | 81 81 16     | Periodic         | Periodic          | MPI + OpenMP                     |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Advection-OMP-xy-ML                    | 80 80 16     | Periodic         | Periodic          | MPI + OpenMP                     |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | multilevel                       |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Advection-xy                           | 81 81 16     | Periodic         | Periodic          | MPI                              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Advection-xy-restart                   | 81 81 16     | Periodic         | Periodic          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | restart                          |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Advection-xy-restart-ML                | 80 80 16     | Periodic         | Periodic          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | restart                          |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | multilevel                       |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| AdvectionU3-xy                         | 81 81 16     | Periodic         | Periodic          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | advection: upstream 3rd order    |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+

Based on :ref:`BioToy<biotoy>`, a small doubly periodic box initialized and forced from NetCDF, which
always includes PnetCDF, Coriolis, non-flat bathymetry, GLS mixing scheme, a nonlinear equation of state,
quadratic bottom stress, bulk fluxes from NetCDF surface forcing, and the Fennel biology model with
carbon, denitrification, and bottom-sediment fluxes:

+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Test                                   | nx ny nz     | xbc              | ybc               | Other                            |
+========================================+==============+==================+===================+==================================+
| BioToy                                 | 4 4 30       | Periodic         | Periodic          | MPI                              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| BioToy-1grid                           | 4 4 30       | Periodic         | Periodic          |                                  |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| BioToy-restart                         | 4 4 30       | Periodic         | Periodic          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | restart                          |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+

Based on :ref:`Boundary Layer<boundarylayer>`, which always includes Coriolis, GLS mixing scheme, non-flat
bathymetry, quadratic bottom stress, bulk fluxes, cloud cover, and evaporation/precipitation with sea
surface height correction:

+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Test                                   | nx ny nz     | xbc              | ybc               | Other                            |
+========================================+==============+==================+===================+==================================+
| BoundaryLayer                          | 39 4 30      | Radiation        | Periodic          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              | / outflow        |                   |                                  |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| BoundaryLayer-OMP                      | 39 4 30      | Radiation        | Periodic          | MPI + OpenMP                     |
|                                        |              |                  |                   |                                  |
|                                        |              | / outflow        |                   |                                  |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| BoundaryLayer-OMP-1grid                | 39 4 30      | Radiation        | Periodic          | OpenMP                           |
|                                        |              |                  |                   |                                  |
|                                        |              | / outflow        |                   |                                  |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| BoundaryLayer-restart                  | 39 4 30      | Radiation        | Periodic          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              | / outflow        |                   | restart                          |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+

Based on :ref:`Channel Test<channeltest>`, which always includes Coriois, GLS mixing scheme, and non-flat bathymetry:

+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Test                                   | nx ny nz     | xbc              | ybc               | Other                            |
+========================================+==============+==================+===================+==================================+
| ChannelTest                            | 20 60 50     | Periodic         | SlipWall          | MPI                              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| ChannelTest-OMP                        | 20 60 50     | Periodic         | SlipWall          | MPI + OpenMP                     |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| ChannelTest-OMP-1grid-xy               | 20 60 50     | Periodic         | SlipWall          | OpenMP                           |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| ChannelTest-xy-restart                 | 20 60 50     | Periodic         | SlipWall          | MPI + OpenMP                     |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | restart                          |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+

Based on :ref:`Dogbone<dogbone>`, which always includes PnetCDF and land-sea masking:

+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| Test                                   | nx ny nz     | xbc              | ybc               | Other                                |
+========================================+==============+==================+===================+======================================+
| Dogbone                                | 42 15 10     | Slipwall         | SlipWall          | MPI                                  |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| Dogbone-1grid                          | 42 15 10     | Slipwall         | SlipWall          |                                      |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| Dogbone-OMP                            | 42 15 10     | Slipwall         | SlipWall          | MPI + OpenMP                         |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| Dogbone-OMP-1grid                      | 42 15 10     | Slipwall         | SlipWall          | OpenMP                               |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| Dogbone-restart                        | 42 15 10     | Slipwall         | Slipwall          | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | restart                              |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+


Based on :ref:`DogboneAnalytic<dogboneanalytic>`, the analytically initialized version of
:ref:`Dogbone<dogbone>`, which always includes land-sea masking, non-flat bathymetry, quadratic bottom
stress, and one refined level at ratio 3. ``MLquad`` refines a fixed lower-left quadrant of the domain,
while ``MLvel`` regrids on the cells where the absolute x-velocity exceeds 0.05:

+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Test                                   | nx ny nz     | xbc              | ybc               | Other                            |
+========================================+==============+==================+===================+==================================+
| DogboneAnalytic_MLquad                 | 42 15 16     | SlipWall         | SlipWall          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | static multilevel                |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DogboneAnalytic_MLquad-1grid           | 42 15 16     | SlipWall         | SlipWall          | static multilevel                |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DogboneAnalytic_MLquad-OMP             | 42 15 16     | SlipWall         | SlipWall          | MPI + OpenMP                     |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | static multilevel                |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DogboneAnalytic_MLquad-OMP-1grid       | 42 15 16     | SlipWall         | SlipWall          | OpenMP                           |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | static multilevel                |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DogboneAnalytic_MLquad-restart         | 42 15 16     | SlipWall         | SlipWall          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | static multilevel                |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | restart                          |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DogboneAnalytic_MLvel                  | 42 15 16     | SlipWall         | SlipWall          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | dynamic multilevel               |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DogboneAnalytic_MLvel-1grid            | 42 15 16     | SlipWall         | SlipWall          | dynamic multilevel               |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DogboneAnalytic_MLvel-OMP              | 42 15 16     | SlipWall         | SlipWall          | MPI + OpenMP                     |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | dynamic multilevel               |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DogboneAnalytic_MLvel-OMP-1grid        | 42 15 16     | SlipWall         | SlipWall          | OpenMP                           |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | dynamic multilevel               |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DogboneAnalytic_MLvel-restart          | 42 15 16     | SlipWall         | SlipWall          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | dynamic multilevel               |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | restart                          |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+

Based on :ref:`Double Gyre<doublegyre>`, which always includes Coriolis:

+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Test                                   | nx ny nz     | xbc              | ybc               | Other                            |
+========================================+==============+==================+===================+==================================+
| DoubleGyre                             | 54 108 4     | SlipWall         | SlipWall          | MPI                              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DoubleGyre-OMP                         | 54 108 4     | SlipWall         | SlipWall          | MPI + OpenMP                     |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DoubleGyre-OMP-1grid-xy                | 54 108 4     | SlipWall         | SlipWall          | OpenMP                           |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DoubleGyre-xy-restart                  | 54 108 4     | SlipWall         | SlipWall          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | restart                          |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+

Based on :ref:`Doubly Periodic<doublyperiodic>`, which always includes Coriolis:

+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Test                                   | nx ny nz     | xbc              | ybc               | Other                            |
+========================================+==============+==================+===================+==================================+
| DoublyPeriodic-1grid-xy                | 41 80 16     | Periodic         | Periodic          |                                  |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DoublyPeriodic-1grid-xy-bathy          | 41 80 16     | Periodic         | Periodic          | non-flat bathyemtry              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DoublyPeriodic-NETCDF-build            | 41 80 16     | N/A              | N/A               | Build w/PnetCDF                  |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DoublyPeriodic-OMP-1grid-xy            | 41 80 16     | Periodic         | Periodic          | OpenMP                           |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DoublyPeriodic-OMP-xy                  | 41 80 16     | Periodic         | Periodic          | MPI + OpenMP                     |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DoublyPeriodic-OMP-xy-bathy            | 41 80 16     | Periodic         | Periodic          | MPI + OpenMP                     |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | non-flat bathymetry              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DoublyPeriodic-xy                      | 41 80 16     | Periodic         | Periodic          | MPI                              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DoublyPeriodic-xy-bathy                | 41 80 16     | Periodic         | Periodic          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | non-flat bathymetry              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DoublyPeriodic-xy-restart              | 41 80 16     | Periodic         | Periodic          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | restart                          |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DoublyPeriodic64-OMP-xy                | 328 320 64   | Periodic         | Periodic          | MPI + OpenMP,                    |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | large problem                    |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DoublyPeriodic64-OMP-xy-bathy          | 328 320 64   | Periodic         | Periodic          | MPI + OpenMP                     |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | large problem                    |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | non-flat bathymetry              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DoublyPeriodic64-xy                    | 328 320 64   | Periodic         | Periodic          | MPI, large problem               |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DoublyPeriodicC4-xy                    | 41 80 16     | Periodic         | Periodic          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | adv.: centered 4th order         |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+

Based on :ref:`Ideal Mini Grid<idealminigrid>`, which always includes Coriolis and PnetCDF. Replace braces in test name with ``IdealMiniGrid``. C-F stands for Chapman-Flather:

+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| Test                                   | nx ny nz     | xbc              | ybc               | Other                                |
+========================================+==============+==================+===================+======================================+
| {}                                     | 10 16 20     | Clamped          | Clamped           | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary salt at boundary                |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}-1grid                               | 10 16 20     | Clamped          | Clamped           | Vary salt at boundary                |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}-CF-Uvel-OMP                         | 10 16 20     | C-F              | C-F               | MPI + OpenMP                         |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary velocity at boundary            |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}-CFO-Salt-OMP                        | 10 16 20     | C-F              | C-F               | MPI + OpenMP                         |
|                                        |              |                  |                   |                                      |
|                                        |              | Radiation        | Radiation         | Vary salt at boundary                |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}-CFO-Temp-OMP                        | 10 16 20     | C-F              | C-F               | MPI + OpenMP                         |
|                                        |              |                  |                   |                                      |
|                                        |              | Radiation        | Radiation         | Vary temperature at boundary         |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}-CFO-Uvel-OMP                        | 10 16 20     | C-F              | C-F               | MPI + OpenMP                         |
|                                        |              |                  |                   |                                      |
|                                        |              | Radiation        | Radiation         | Vary velocity at boundary            |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}-EWWall-OMP                          | 10 16 20     | SlipWall         | Clamped           | MPI + OpenMP                         |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary salt at boundary                |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}-EWWall-restart                      | 10 16 20     | SlipWall         | Clamped           | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary salt at boundary                |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | restart                              |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}-NSWall-OMP                          | 10 16 20     | Clamped          | SlipWall          | MPI + OpenMP                         |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary salt at boundary                |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}-NSWall-restart                      | 10 16 20     | Clamped          | Slipwall          | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary salt at boundary                |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | restart                              |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}-OMP                                 | 10 16 20     | Clamped          | Clamped           | MPI + OpenMP                         |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary salt at boundary                |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}-OMP-1grid                           | 10 16 20     | Clamped          | Clamped           | OpenMP                               |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary salt at boundary                |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}-Temp                                | 10 16 20     | Clamped          | Clamped           | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary temperature at boundary         |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}-Uvel                                | 10 16 20     | Clamped          | Clamped           | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary velocity at boundary            |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}-Uvel-1grid                          | 10 16 20     | Clamped          | Clamped           | Vary velocity at boundary            |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}-Uvel-EWWall-OMP                     | 10 16 20     | SlipWall         | Clamped           | MPI + OpenMP                         |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary velocity at boundary            |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}-Uvel-NSWall-OMP                     | 10 16 20     | Clamped          | SlipWall          | MPI + OpenMP                         |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary velocity at boundary            |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}-Uvel-OMP                            | 10 16 20     | Clamped          | Clamped           | MPI + OpenMP                         |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary velocity at boundary            |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}-Uvel-OMP-1grid                      | 10 16 20     | Clamped          | Clamped           | OpenMP                               |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary velocity at boundary            |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}-restart                             | 10 16 20     | Clamped          | Clamped           | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary salt at boundary                |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | restart                              |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}_synth_hires                         | 10 10 4      | Outflow          | Outflow           | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | multilevel (ratio 3)                 |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | high-resolution grid and IC          |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | no Coriolis                          |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}_synth_hires-1grid                   | 10 10 4      | Outflow          | Outflow           | multilevel (ratio 3)                 |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | high-resolution grid and IC          |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | no Coriolis                          |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}_synth_hires-restart                 | 10 10 4      | Outflow          | Outflow           | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | multilevel (ratio 3)                 |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | high-resolution grid and IC          |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | no Coriolis                          |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | restart                              |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Mask                                 | 10 16 20     | Clamped          | Clamped           | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary salt at boundary                |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | land-sea masking                     |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Mask-1grid                           | 10 16 20     | Clamped          | Clamped           | Coriolis                             |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | land-sea masking                     |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Mask-CF-Uvel-OMP                     | 10 16 20     | C-F              | C-F               | MPI + OpenMP                         |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary velocity at boundary            |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | land-sea masking                     |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Mask-CFO-Salt-OMP                    | 10 16 20     | C-F              | C-F               | MPI + OpenMP                         |
|                                        |              |                  |                   |                                      |
|                                        |              | Radiation        | Radiaion          | Vary salt at boundary                |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | land-sea masking                     |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Mask-CFO-Temp-OMP                    | 10 16 20     | C-F              | C-F               | MPI + OpenMP                         |
|                                        |              |                  |                   |                                      |
|                                        |              | Radiation        | Radiaion          | Vary temperature at boundary         |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | land-sea masking                     |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Mask-CFO-Uvel-OMP                    | 10 16 20     | C-F              | C-F               | MPI + OpenMP                         |
|                                        |              |                  |                   |                                      |
|                                        |              | Radiation        | Radiaion          | Vary velocity at boundary            |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | land-sea masking                     |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Mask-EWWall-OMP                      | 10 16 20     | SlipWall         | Clamped           | MPI + OpenMP                         |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary salt at boundary                |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | land-sea masking                     |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Mask-NSWall-OMP                      | 10 16 20     | Clamped          | SlipWall          | MPI + OpenMP                         |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary salt at boundary                |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | land-sea masking                     |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Mask-OMP                             | 10 16 20     | Clamped          | Clamped           | MPI + OpenMP                         |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary salt at boundary                |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | land-sea masking                     |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Mask-OMP-1grid                       | 10 16 20     | Clamped          | Clamped           | OpenMP                               |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary salt at boundary                |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | land-sea masking                     |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Mask-Temp                            | 10 16 20     | Clamped          | Clamped           | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary temperature at boundary         |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | land-sea masking                     |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Mask-Uvel                            | 10 16 20     | Clamped          | Clamped           | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary velocity at boundary            |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | land-sea masking                     |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Mask-Uvel-1grid                      | 10 16 20     | Clamped          | Clamped           | Vary velocity at boundary            |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | land-sea masking                     |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Mask-Uvel-EWWall-OMP                 | 10 16 20     | SlipWall         | Clamped           | MPI + OpenMP                         |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary velocity at boundary            |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | land-sea masking                     |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Mask-Uvel-NSWall-OMP                 | 10 16 20     | Clamped          | SlipWall          | MPI + OpenMP                         |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary velocity at boundary            |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | land-sea masking                     |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Mask-Uvel-OMP                        | 10 16 20     | Clamped          | Clamped           | MPI + OpenMP                         |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary velocity at boundary            |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | land-sea masking                     |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Mask-Uvel-OMP-1grid                  | 10 16 20     | Clamped          | Clamped           | OpenMP                               |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary velocity at boundary            |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | land-sea masking                     |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Mask-restart                         | 10 16 20     | Clamped          | Clamped           | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary salt at boundary                |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | land-sea masking                     |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | restart                              |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Wind                                 | 10 16 20     | Clamped          | Clamped           | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary salt at boundary                |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Surface wind and bulk fluxes         |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Wind-1grid                           | 10 16 20     | Clamped          | Clamped           | Vary salt at boundary                |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Surface wind and bulk fluxes         |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Wind-OMP                             | 10 16 20     | Clamped          | Clamped           | MPI + OpenMP                         |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary salt at boundary                |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Surface wind and bulk fluxes         |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Wind-OMP-1grid                       | 10 16 20     | Clamped          | Clamped           | OpenMP                               |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary salt at boundary                |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Surface wind and bulk fluxes         |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Wind-restart                         | 10 16 20     | Clamped          | Clamped           | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary salt at boundary                |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | restart                              |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Surface wind and bulk fluxes         |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Clim                                 | 10 16 20     | C-F              | C-F               | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              | Radiation        | Radiation         | Vary salt, temp, v at boundary       |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Climatology nudging                  |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Clim-OMP                             | 10 16 20     | C-F              | C-F               | MPI + OpenMP                         |
|                                        |              |                  |                   |                                      |
|                                        |              | Radiation        | Radiation         | Vary salt, temp, v at boundary       |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Climatology nudging                  |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Clim-OMP-1grid                       | 10 16 20     | C-F              | C-F               | OpenMP                               |
|                                        |              |                  |                   |                                      |
|                                        |              | Radiation        | Radiation         | Vary salt, temp, v at boundary       |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Climatology nudging                  |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Clim-restart                         | 10 16 20     | C-F              | C-F               | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              | Radiation        | Radiation         | Vary salt, temp, v at boundary       |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | restart                              |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Climatology nudging                  |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+

Based on :ref:`Ideal River Grid<idealrivgrid>`, which always includes Coriolis, PnetCDF, and rivers:

+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| Test                                   | nx ny nz     | xbc              | ybc               | Other                                |
+========================================+==============+==================+===================+======================================+
| IdealRivGrid                           | 10 16 20     | Clamped          | Clamped           | MPI                                  |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| IdealRivGrid-1grid                     | 10 16 20     | Clamped          | Clamped           |                                      |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| IdealRivGrid-OMP                       | 10 16 20     | Clamped          | Clamped           | MPI + OpenMP                         |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| IdealRivGrid-OMP-1grid                 | 10 16 20     | Clamped          | Clamped           | OpenMP                               |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| IdealRivGrid-restart                   | 10 16 20     | Clamped          | Clamped           | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | restart                              |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+

Based on :ref:`Particle Advection Flat<particleadvectionflat>`, which always include MPI, Coriolis, and tracer particles:

+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Test                                   | nx ny nz     | xbc              | ybc               | Other                            |
+========================================+==============+==================+===================+==================================+
| ParticleAdvectionFlat                  | 41 80 16     | Periodic         | Periodic          |                                  |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| ParticleAdvectionFlat-restart          | 41 80 16     | Periodic         | Periodic          | restart                          |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+

Based on :ref:`Seamount<seamount-desc>`, which always includes Coriolis and non-flat bathymetry:

+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Test                                   | nx ny nz     | xbc              | ybc               | Other                            |
+========================================+==============+==================+===================+==================================+
| Seamount-1grid-xy                      | 49 48 13     | Periodic         | Periodic          |                                  |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Seamount-OMP-xy                        | 49 48 13     | Periodic         | Periodic          | MPI + OpenMP                     |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Seamount-xy                            | 49 48 13     | Periodic         | Periodic          | MPI                              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Seamount-xy-restart                    | 49 48 13     | Periodic         | Periodic          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | restart                          |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Seamount64-OMP-xy                      | 320 320 64   | Periodic         | Periodic          | MPI + OpenMP, large problem      |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+

Based on :ref:`Upwelling<upwelling-desc>`, which always includes Coriolis and non-flat bathymetry:

+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Test                                   | nx ny nz     | xbc              | ybc               | Other                            |
+========================================+==============+==================+===================+==================================+
| Upwelling                              | 41 80 16     | Periodic         | SlipWall          | MPI                              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling-1grid                        | 41 80 16     | Periodic         | SlipWall          | Coriolis                         |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling-OMP                          | 41 80 16     | Periodic         | SlipWall          | MPI + OpenMP                     |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling-logDrag-OMP                  | 41 80 16     | Periodic         | SlipWall          | MPI + OpenMP                     |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | logarithmic bottom stress        |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling-nonlinEOS-OMP                | 41 80 16     | Periodic         | SlipWall          | MPI + OpenMP                     |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | nonlinear equation of state      |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling-quadDrag-OMP                 | 41 80 16     | Periodic         | SlipWall          | MPI + OpenMP                     |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | quadratic bottom stress          |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling-OMP-1grid                    | 41 80 16     | Periodic         | SlipWall          | OpenMP                           |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling-restart                      | 41 80 16     | Periodic         | SlipWall          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | restart                          |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling-x                            | 41 80 16     | SlipWall         | Periodic          | MPI                              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling-x-1grid                      | 41 80 16     | SlipWall         | Periodic          |                                  |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling-x-OMP                        | 41 80 16     | SlipWall         | Periodic          | MPI + OpenMP                     |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling64-OMP                        | 328 320 64   | SlipWall         | Periodic          | MPI + OpenMP, large problem      |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| UpwellingC4                            | 41 80 16     | Periodic         | SlipWall          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | advection: centered 4th order    |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling_GeopotentialMixing           | 41 80 16     | Periodic         | SlipWall          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | constant horizontal mixing       |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | harmonic tracer diffusion        |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | rotated along geopotential       |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling_GLS                          | 41 80 16     | Periodic         | SlipWall          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | GLS mixing scheme                |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling_GLS-restart                  | 41 80 16     | Periodic         | SlipWall          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | GLS mixing scheme                |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | restart                          |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling_GLS_Canuto_A                 | 41 80 16     | Periodic         | SlipWall          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | GLS mixing scheme                |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | Canuto A stability               |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling_GLS_Canuto_B                 | 41 80 16     | Periodic         | SlipWall          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | GLS mixing scheme                |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | Canuto B stability               |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+

Based on :ref:`Upwelling<upwelling-desc>` with the Fennel biology model, which always includes Coriolis,
non-flat bathymetry, bulk fluxes, and analytically initialized Fennel biology. The tests without ``hires``
additionally enable carbon, oxygen, river DON, non-conservative alkalinity, a dated atmospheric pCO2, and
Wanninkhof (2014) gas transfer. The ``hires`` tests instead add one refined level at ratio 3, with the
level-0 bathymetry averaged down from it:

+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Test                                   | nx ny nz     | xbc              | ybc               | Other                            |
+========================================+==============+==================+===================+==================================+
| Upwelling-Fennel                       | 41 80 16     | Periodic         | SlipWall          | MPI                              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling-Fennel-1grid                 | 41 80 16     | Periodic         | SlipWall          |                                  |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling-Fennel-hires                 | 41 80 16     | Periodic         | SlipWall          | MPI                              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling-Fennel-hires-1grid           | 41 80 16     | Periodic         | SlipWall          |                                  |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling-fennel-hires-restart         | 41 80 16     | Periodic         | SlipWall          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | restart                          |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling-fennel-restart               | 41 80 16     | Periodic         | SlipWall          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | restart                          |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+

Nightly Regression Tests on GPU
-------------------------------

And the following are currently tested nighly on GPU. All are compiled and run with CUDA.

Based on :ref:`Advection<advection>`:

+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Test                                   | nx ny nz     | xbc              | ybc               | Other                            |
+========================================+==============+==================+===================+==================================+
| Advection-1grid-xy                     | 81 81 16     | Periodic         | Periodic          |                                  |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Advection-1grid-xy-ML                  | 80 80 16     | Periodic         | Periodic          | multilevel                       |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Advection-xy                           | 81 81 16     | Periodic         | Periodic          | MPI                              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Advection-xy-ML                        | 80 80 16     | Periodic         | Periodic          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | multilevel                       |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Advection-xy-restart                   | 81 81 16     | Periodic         | Periodic          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | restart                          |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Advection64-xy                         | 328 320 64   | Periodic         | Periodic          | MPI, large problem               |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | restart                          |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| AdvectionU3-xy                         | 81 81 16     | Periodic         | Periodic          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | advection: upstream 3rd order    |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+

Based on :ref:`BioToy<biotoy>`, a small doubly periodic box initialized and forced from NetCDF, which
always includes PnetCDF, Coriolis, non-flat bathymetry, GLS mixing scheme, a nonlinear equation of state,
quadratic bottom stress, bulk fluxes from NetCDF surface forcing, and the Fennel biology model with
carbon, denitrification, and bottom-sediment fluxes:

+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Test                                   | nx ny nz     | xbc              | ybc               | Other                            |
+========================================+==============+==================+===================+==================================+
| BioToy                                 | 4 4 30       | Periodic         | Periodic          | MPI                              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| BioToy-1grid                           | 4 4 30       | Periodic         | Periodic          |                                  |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| BioToy-restart                         | 4 4 30       | Periodic         | Periodic          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | restart                          |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+

Based on :ref:`Boundary Layer<boundarylayer>`, which always includes Coriolis, GLS mixing scheme, non-flat
bathymetry, quadratic bottom stress, bulk fluxes, cloud cover, and evaporation/precipitation with sea
surface height correction. This problem is sensitive enough on GPU that its comparison
tolerance is loosened to ``2e-5``:

+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Test                                   | nx ny nz     | xbc              | ybc               | Other                            |
+========================================+==============+==================+===================+==================================+
| BoundaryLayer-1grid-xy                 | 39 4 30      | Radiation        | Periodic          |                                  |
|                                        |              |                  |                   |                                  |
|                                        |              | / outflow        |                   |                                  |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| BoundaryLayer-xy                       | 39 4 30      | Radiation        | Periodic          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              | / outflow        |                   |                                  |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| BoundaryLayer-xy-restart               | 39 4 30      | Radiation        | Periodic          | restart                          |
|                                        |              |                  |                   |                                  |
|                                        |              | / outflow        |                   |                                  |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+

Based on :ref:`Channel Test<channeltest>`, which always includes Coriolis, GLS mixing scheme, and non-flat bathymetry:

+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Test                                   | nx ny nz     | xbc              | ybc               | Other                            |
+========================================+==============+==================+===================+==================================+
| ChannelTest-1grid-xy                   | 20 60 50     | Periodic         | SlipWall          |                                  |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| ChannelTest-xy                         | 20 60 50     | Periodic         | SlipWall          | MPI                              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| ChannelTest-xy-restart                 | 20 60 50     | Periodic         | SlipWall          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | restart                          |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+

Based on :ref:`Dogbone<dogbone>`, which always includes PnetCDF and land-sea masking:

+----------------------------------------+--------------+------------------+-------------------+-------------+
| Test                                   | nx ny nz     | xbc              | ybc               | Other       |
+========================================+==============+==================+===================+=============+
| Dogbone                                | 42 15 10     | Slipwall         | Slipwall          | MPI         |
+----------------------------------------+--------------+------------------+-------------------+-------------+
| Dogbone-1grid                          | 42 15 10     | Slipwall         | Slipwall          |             |
+----------------------------------------+--------------+------------------+-------------------+-------------+
| Dogbone-restart                        | 42 15 10     | Slipwall         | Slipwall          | MPI         |
|                                        |              |                  |                   |             |
|                                        |              |                  |                   | restart     |
+----------------------------------------+--------------+------------------+-------------------+-------------+


Based on :ref:`DogboneAnalytic<dogboneanalytic>`, the analytically initialized version of
:ref:`Dogbone<dogbone>`, which always includes land-sea masking, non-flat bathymetry, quadratic bottom
stress, and one refined level at ratio 3. ``MLquad`` refines a fixed lower-left quadrant of the domain,
while ``MLvel`` regrids on the cells where the absolute x-velocity exceeds 0.05:

+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Test                                   | nx ny nz     | xbc              | ybc               | Other                            |
+========================================+==============+==================+===================+==================================+
| DogboneAnalytic_MLquad                 | 42 15 16     | SlipWall         | SlipWall          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | static multilevel                |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DogboneAnalytic_MLquad-1grid           | 42 15 16     | SlipWall         | SlipWall          | static multilevel                |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DogboneAnalytic_MLquad-restart         | 42 15 16     | SlipWall         | SlipWall          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | static multilevel                |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | restart                          |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DogboneAnalytic_MLvel                  | 42 15 16     | SlipWall         | SlipWall          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | dynamic multilevel               |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DogboneAnalytic_MLvel-1grid            | 42 15 16     | SlipWall         | SlipWall          | dynamic multilevel               |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DogboneAnalytic_MLvel-restart          | 42 15 16     | SlipWall         | SlipWall          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | dynamic multilevel               |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | restart                          |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+

Based on :ref:`Double Gyre<doublegyre>`, which always includes Coriolis:

+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Test                                   | nx ny nz     | xbc              | ybc               | Other                            |
+========================================+==============+==================+===================+==================================+
| DoubleGyre-1grid-xy                    | 54 108 4     | SlipWall         | SlipWall          |                                  |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DoubleGyre-xy                          | 54 108 4     | SlipWall         | SlipWall          | MPI                              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DoubleGyre-xy-restart                  | 54 108 4     | SlipWall         | SlipWall          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | restart                          |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+

Based on :ref:`Doubly Periodic<doublyperiodic>`, which always includes Coriolis:

+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Test                                   | nx ny nz     | xbc              | ybc               | Other                            |
+========================================+==============+==================+===================+==================================+
| DoublyPeriodic-1grid-xy                | 41 80 16     | Periodic         | Periodic          |                                  |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DoublyPeriodic-xy                      | 41 80 16     | Periodic         | Periodic          | MPI                              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DoublyPeriodic-xy-bathy                | 41 80 16     | Periodic         | Periodic          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | non-flat bathymetry              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DoublyPeriodic-xy-restart              | 41 80 16     | Periodic         | Periodic          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | restart                          |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DoublyPeriodic64-xy                    | 328 320 64   | Periodic         | Periodic          | MPI, large problem               |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DoublyPeriodic64-xy-bathy              | 328 320 64   | Periodic         | Periodic          | MPI, large problem               |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | non-flat bathymetry              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| DoublyPeriodicC4-xy                    | 41 80 16     | Periodic         | Periodic          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | adv.: centered 4th order         |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+

Based on :ref:`Ideal Mini Grid<idealminigrid>`, which always includes Coriolis and PnetCDF. Replace braces in test name with ``IdealMiniGrid``. C-F stands for Chapman-Flather:

+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| Test                                   | nx ny nz     | xbc              | ybc               | Other                                |
+========================================+==============+==================+===================+======================================+
| {}                                     | 10 16 20     | Clamped          | Clamped           | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary salt at boundary                |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}-1grid                               | 10 16 20     | Clamped          | Clamped           | Vary salt at boundary                |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}-CF-Uvel                             | 10 16 20     | C-F              | C-F               | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary velocity at boundary            |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}-CFO-Salt                            | 10 16 20     | C-F              | C-F               | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              | Radiation        | Radiation         | Vary salt at boundary                |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}-CFO-Temp                            | 10 16 20     | C-F              | C-F               | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              | Radiation        | Radiation         | Vary temperature at boundary         |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}-CFO-Uvel                            | 10 16 20     | C-F              | C-F               | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              | Radiation        | Radiation         | Vary velocity at boundary            |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}-EWWall                              | 10 16 20     | SlipWall         | Clamped           | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary salt at boundary                |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}-NSWall                              | 10 16 20     | Clamped          | SlipWall          | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary salt at boundary                |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}-Temp                                | 10 16 20     | Clamped          | Clamped           | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary temperature at boundary         |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}-Uvel                                | 10 16 20     | Clamped          | Clamped           | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary velocity at boundary            |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}-Uvel-1grid                          | 10 16 20     | Clamped          | Clamped           | Vary velocity at boundary            |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}-Uvel-EWWall                         | 10 16 20     | SlipWall         | Clamped           | Vary velocity at boundary            |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}-Uvel-NSWall                         | 10 16 20     | Clamped          | SlipWall          | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary velocity at boundary            |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}-restart                             | 10 16 20     | Clamped          | Clamped           | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary salt at boundary                |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | restart                              |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}_synth_hires                         | 10 10 4      | Outflow          | Outflow           | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | multilevel (ratio 3)                 |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | high-resolution grid and IC          |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | no Coriolis                          |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}_synth_hires-1grid                   | 10 10 4      | Outflow          | Outflow           | multilevel (ratio 3)                 |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | high-resolution grid and IC          |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | no Coriolis                          |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}_synth_hires-restart                 | 10 10 4      | Outflow          | Outflow           | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | multilevel (ratio 3)                 |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | high-resolution grid and IC          |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | no Coriolis                          |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | restart                              |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Mask                                 | 10 16 20     | Clamped          | Clamped           | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary salt at boundary                |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | land-sea masking                     |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Mask-1grid                           | 10 16 20     | Clamped          | Clamped           | Vary salt at boundary                |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | land-sea masking                     |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Mask-CF-Uvel                         | 10 16 20     | C-F              | C-F               | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary velocity at boundary            |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | land-sea masking                     |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Mask-CFO-Salt                        | 10 16 20     | C-F              | C-F               | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              | Radiation        | Radiaion          | Vary salt at boundary                |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | land-sea masking                     |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Mask-CFO-Temp                        | 10 16 20     | C-F              | C-F               | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              | Radiation        | Radiaion          | Vary temperature at boundary         |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | land-sea masking                     |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Mask-CFO-Uvel                        | 10 16 20     | C-F              | C-F               | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              | Radiation        | Radiaion          | Vary velocity at boundary            |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | land-sea masking                     |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Mask-EWWall                          | 10 16 20     | SlipWall         | Clamped           | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary salt at boundary                |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | land-sea masking                     |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Mask-EWWall-restart                  | 10 16 20     | SlipWall         | Clamped           | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary salt at boundary                |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | land-sea masking                     |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | restart                              |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Mask-NSWall                          | 10 16 20     | Clamped          | Slipwall          | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary salt at boundary                |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | land-sea masking                     |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Mask-NSWall-restart                  | 10 16 20     | Clamped          | Slipwall          | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary salt at boundary                |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | land-sea masking                     |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | restart                              |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Mask-Temp                            | 10 16 20     | Clamped          | Clamped           | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary temperature at boundary         |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | land-sea masking                     |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Mask-Uvel                            | 10 16 20     | Clamped          | Clamped           | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary velocity at boundary            |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | land-sea masking                     |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Mask-Uvel-1grid                      | 10 16 20     | Clamped          | Clamped           | Vary velocity at boundary            |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | land-sea masking                     |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Mask-Uvel-EWWall                     | 10 16 20     | SlipWall         | Clamped           | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary salt at boundary                |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | land-sea masking                     |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Mask-Uvel-NSWall                     | 10 16 20     | Clamped          | Periodic          | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary salt at boundary                |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | land-sea masking                     |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Mask-restart                         | 10 16 20     | Clamped          | Clamped           | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary salt at boundary                |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | land-sea masking                     |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | restart                              |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Wind                                 | 10 16 20     | Clamped          | Clamped           | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary salt at boundary                |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Surface wind and bulk fluxes         |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Wind-1grid                           | 10 16 20     | Clamped          | Clamped           | Vary salt at boundary                |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Surface wind and bulk fluxes         |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Wind-restart                         | 10 16 20     | Clamped          | Clamped           | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Vary salt at boundary                |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | restart                              |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Surface wind and bulk fluxes         |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Clim                                 | 10 16 20     | C-F              | C-F               | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              | Radiation        | Radiation         | Vary salt, temp, v at boundary       |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Climatology nudging                  |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Clim-1grid                           | 10 16 20     | C-F              | C-F               | Vary salt, temp, v at boundary       |
|                                        |              |                  |                   |                                      |
|                                        |              | Radiation        | Radiation         | Climatology nudging                  |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| {}Clim-restart                         | 10 16 20     | C-F              | C-F               | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              | Radiation        | Radiation         | Vary salt, temp, v at boundary       |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | restart                              |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | Climatology nudging                  |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+

Based on :ref:`Ideal River Grid<idealrivgrid>`, which always includes Coriolis, PnetCDF, and rivers:

+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| Test                                   | nx ny nz     | xbc              | ybc               | Other                                |
+========================================+==============+==================+===================+======================================+
| IdealRivGrid                           | 10 16 20     | Clamped          | Clamped           | MPI                                  |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| IdealRivGrid-1grid                     | 10 16 20     | Clamped          | Clamped           |                                      |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+
| IdealRivGrid-restart                   | 10 16 20     | Clamped          | Clamped           | MPI                                  |
|                                        |              |                  |                   |                                      |
|                                        |              |                  |                   | restart                              |
+----------------------------------------+--------------+------------------+-------------------+--------------------------------------+

Based on :ref:`Particle Advection Flat<particleadvectionflat>`, which always includes MPI, Coriolis, and tracer particles:

+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Test                                   | nx ny nz     | xbc              | ybc               | Other                            |
+========================================+==============+==================+===================+==================================+
| ParticleAdvectionFlat                  | 41 80 16     | Periodic         | Periodic          |                                  |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| ParticleAdvectionFlat-restart          | 41 80 16     | Periodic         | Periodic          | restart                          |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+

Based on :ref:`Seamount<seamount-desc>`, which always includes Coriolis and non-flat bathymetry:

+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Test                                   | nx ny nz     | xbc              | ybc               | Other                            |
+========================================+==============+==================+===================+==================================+
| Seamount-1grid-xy                      | 49 48 13     | Periodic         | Periodic          |                                  |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Seamount-xy                            | 49 48 13     | Periodic         | Periodic          | MPI                              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Seamount64-xy                          | 320 320 64   | Periodic         | Periodic          | MPI, large problem               |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+

Based on :ref:`Upwelling<upwelling-desc>`, which always includes Coriolis and non-flat bathymetry:

+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Test                                   | nx ny nz     | xbc              | ybc               | Other                            |
+========================================+==============+==================+===================+==================================+
| Upwelling                              | 41 80 16     | Periodic         | SlipWall          | MPI                              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling-1grid                        | 41 80 16     | Periodic         | SlipWall          |                                  |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling-logDrag                      | 41 80 16     | Periodic         | SlipWall          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | logarithmic bottom stress        |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling-nonlinEOS                    | 41 80 16     | Periodic         | SlipWall          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | nonlinear equation of state      |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling-quadDrag                     | 41 80 16     | Periodic         | SlipWall          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | quadratic bottom stress          |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling-x                            | 41 80 16     | SlipWall         | Periodic          | MPI                              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling-x-1grid                      | 41 80 16     | SlipWall         | Periodic          |                                  |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling64                            | 328 320 64   | SlipWall         | Periodic          | MPI, large problem               |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| UpwellingC4                            | 41 80 16     | Periodic         | SlipWall          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | advection: centered 4th order    |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling-GeopotentialMixing           | 41 80 16     | Periodic         | SlipWall          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | constant horizontal mixing       |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | harmonic tracer diffusion        |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | rotated along geopotential       |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling_GLS                          | 41 80 16     | Periodic         | SlipWall          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | GLS mixing scheme                |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling_GLS-restart                  | 41 80 16     | Periodic         | SlipWall          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | GLS mixing scheme                |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | restart                          |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling_GLS_Canuto_A                 | 41 80 16     | Periodic         | SlipWall          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | GLS mixing scheme                |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | Canuto A stability               |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling_GLS_Canuto_B                 | 41 80 16     | Periodic         | SlipWall          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | GLS mixing scheme                |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | Canuto B stability               |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+

Based on :ref:`Upwelling<upwelling-desc>` with the Fennel biology model, which always includes Coriolis,
non-flat bathymetry, bulk fluxes, and analytically initialized Fennel biology. The tests without ``hires``
additionally enable carbon, oxygen, river DON, non-conservative alkalinity, a dated atmospheric pCO2, and
Wanninkhof (2014) gas transfer. The ``hires`` tests instead add one refined level at ratio 3, with the
level-0 bathymetry averaged down from it:

+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Test                                   | nx ny nz     | xbc              | ybc               | Other                            |
+========================================+==============+==================+===================+==================================+
| Upwelling-Fennel                       | 41 80 16     | Periodic         | SlipWall          | MPI                              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling-Fennel-1grid                 | 41 80 16     | Periodic         | SlipWall          |                                  |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling-Fennel-hires                 | 41 80 16     | Periodic         | SlipWall          | MPI                              |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling-Fennel-hires-1grid           | 41 80 16     | Periodic         | SlipWall          |                                  |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling-Fennel-hires-restart         | 41 80 16     | Periodic         | SlipWall          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | restart                          |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
| Upwelling-Fennel-restart               | 41 80 16     | Periodic         | SlipWall          | MPI                              |
|                                        |              |                  |                   |                                  |
|                                        |              |                  |                   | restart                          |
+----------------------------------------+--------------+------------------+-------------------+----------------------------------+
