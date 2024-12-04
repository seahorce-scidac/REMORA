
 .. role:: cpp(code)
    :language: c++

 .. _MeshRefinement:

Mesh Refinement
===============

REMORA allows both static and dynamic mesh refinement, as well as the choice of one-way or two-way coupling.

Note that any tagged region will be covered by one or more boxes.  The user may
specify the refinement criteria and/or region to be covered, but not the decomposition of the region into
individual grids.

See the `Gridding`_ section of the AMReX documentation for details of how individual grids are created.

.. _`Gridding`: https://amrex-codes.github.io/amrex/docs_html/ManagingGridHierarchy_Chapter.html

Static Mesh Refinement
----------------------

For static refinement, we control the placement of grids by specifying
the low and high extents (in physical space) of each box in the lateral
directions. REMORA enforces that all refinement spans the entire vertical direction.

The following example demonstrates how to tag regions for static refinement.
In this first example, all cells in the region :math:`[0.15,0.25,\texttt{prob_lo_z}] \times [0.35,0.45,\texttt{prob_hi_z}]`
and in the region :math:`[0.65,0.75,\texttt{prob_lo_z}]\times[0.85,0.95,\texttt{prob_hi_z}]` are tagged for
one level of refinement, where prob_lo_z and prob_hi_z are the vertical extents of the domain:

::

          amr.max_level = 1
          amr.ref_ratio = 2

          remora.refinement_indicators = box1 box2

          remora.box1.in_box_lo = .15 .25
          remora.box1.in_box_hi = .35 .45

          remora.box2.in_box_lo = .65 .75
          remora.box2.in_box_hi = .85 .95

In the example below, we refine the region :math:`[0.15,0.25,\texttt{prob_lo_z}]\times [0.35,0.45,\texttt{prob_hi_z}]`
by two levels of factor 3 refinement. In this case, the refined region at level 1 will
be sufficient to enclose the refined region at level 2.

::

          amr.max_level = 2
          amr.ref_ratio = 3 3

          remora.refinement_indicators = box1

          remora.box1.in_box_lo = .15 .25
          remora.box1.in_box_hi = .35 .45

And in this final example, the region :math:`[0.15,0.25,\texttt{prob_lo_z}]\times[0.35,0.45,\texttt{prob_hi_z}]`
will be refined by two levels of factor 3, but the larger region, :math:`[0.05,0.05,\texttt{prob_lo_z}]\times [0.75,0.75,\texttt{prob_hi_z}]``
will be refined by a single factor 3 refinement.

::

          amr.max_level = 2
          amr.ref_ratio = 3 3

          remora.refinement_indicators = box1 box2

          remora.box1.in_box_lo = .15 .25
          remora.box1.in_box_hi = .35 .45

          remora.box2.in_box_lo = .05 .05
          remora.box2.in_box_hi = .75 .75
          remora.box2.max_level = 1


Dynamic Mesh Refinement
-----------------------

Dynamically created tagging functions are based on runtime data specified in the inputs file.
These dynamically generated functions test on either state variables or derived variables
defined in REMORA_derive.cpp and included in the derive_list in Setup.cpp.

Available tests include

-  “greater\_than”: :math:`\text{field} \geq \text{threshold}`

-  “less\_than”: :math:`\text{field} \leq \text{threshold}`

-  “adjacent\_difference\_greater”: :math:`\text{max}( | \text{difference between any nearest-neighbor cell} | ) \geq \text{threshold}`

The example below adds two user-named criteria:

- ``hi_temp``: cells with density greater than 10 on level 0, and greater than 20 on level 1 and higher;
- ``lo_vort``: cells with relative vorticity less than 0 that are inside the region :math:`[0.25,0.25,\texttt{prob_lo_z}]\times[0.75,0.75,\texttt{prob_hi_z}]`;
- ``scalardiff``: cells having a difference in the scalar of 0.01 or more from that of any immediate neighbor.

The first will trigger up to AMR level 3 and the second to level 2.
The second will be active only when the problem time is between 100 and 300 seconds.

Note that ``temp`` and ``scalar`` are the names of state variables and ``vorticity`` is a derived variable.
Valid field options for refinement are: ``scalar``, ``temp``, ``salt``, ``x_velocity``, ``y_velocity``, ``z_velocity``,
and ``vorticity``.

::

          remora.refinement_indicators = hi_temp scalardiff

          remora.hi_temp.max_level = 3
          remora.hi_temp.value_greater = 10. 20.
          remora.hi_temp.field_name = temp

          remora.scalardiff.max_level = 2
          remora.scalardiff.adjacent_difference_greater = 0.01
          remora.scalardiff.field_name = scalar
          remora.scalardiff.start_time = 100
          remora.scalardiff.end_time = 300

          remora.lo_vort.max_level = 1
          remora.lo_vort.value_less = 0
          remora.lo_vort.field_name = vorticity
          remora.lo_vort.in_box_lo = .25 .25
          remora.lo_vort.in_box_hi = .75 .75

Coupling Types
--------------

REMORA supports one-way and two-way coupling between levels; this is a run-time input

::

      remora.coupling_type = "OneWay" or "TwoWay"

By one-way coupling, we mean that between each pair of refinement levels,
the coarse level communicates data to the fine level to serve as boundary conditions
for the time advance of the fine solution. For cell-centered quantities,
and face-baced normal momenta on the coarse-fine interface, the coarse data is conservatively
interpolated to the fine level.

The interpolated data is utilized to specify ghost cell data (outside of the valid fine region).

By two-way coupling, we mean that in additional to interpolating data from the coarser level
to supply boundary conditions for the fine regions,
the fine level also communicates data back to the coarse level in two ways:

- The fine cell-centered data are conservatively averaged onto the coarse mesh covered by fine mesh.

- The fine momenta are conservatively averaged onto the coarse faces covered by fine mesh.

- A "reflux" operation is performed for all cell-centered data; this updates values on the coarser level outside of regions covered by the finer level.

Advected quantities which are advanced in conservation form will lose conservation with one-way coupling.
Two-way coupling ensures conservation of the advective contribution to all scalar updates but
does not account for loss of conservation due to diffusive or source terms.

.. _sec:fillghost:

Filling Ghost Values
--------------------

REMORA uses an operation called ``FillPatch`` to fill the ghost cells/faces for each grid of data.
The data is filled outside the valid region with a combination of three operations: interpolation
from coarser level, copy from same level, and enforcement of physical boundary conditions.

Interpolation from Coarser level
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Interpolation is controlled by which interpolater we choose to use. The default is
conservative interpolation for cell-centered quantities, and analogous for faces.
These options are currently hard-coded in REMORA.
The paradigm is that fine faces on a coarse-fine boundary are filled as Dirichlet
boundary conditions from the coarser level; all faces outside the valid region are
similarly filled, while fine faces inside the valid region are not over-written.

Copy from other grids at same level (includes periodic boundaries)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is part of the ``FillPatch`` operation, but can also be applied independently,
e.g. by the call

::

    mf.FillBoundary(geom[lev].periodicity());

would fill all the ghost cells/faces of the grids in MultiFab ``mf``, including those
that occur at periodic boundaries.

In the ``FillPatch`` operation, ``FillBoundary`` always overrides any interpolated values, i.e. if
there is fine data available (except at coarse-fine boundary) we always use it.

