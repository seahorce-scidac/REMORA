
 .. role:: cpp(code)
    :language: c++

 .. _MeshRefinement:

Mesh Refinement
===============

REMORA allows both static and dynamic mesh refinement, as well as the choice of one-way or two-way coupling.

Note that any tagged region will be covered by one or more boxes.  The user may
specify the refinement criteria and/or region to be covered, but not the decomposition of the region into
individual grids. REMORA enforces that all refinement spans the entire vertical direction. Field-based
tagging criteria are evaluated only where the field means something, so the land-sea boundary does not
drive refinement on its own account (except for a criterion keyed on ``mask`` itself, which exists to
find the coast); a region named explicitly, by a box, is refined in full whether it
is land or water. See `Masked Regions and Tagging`_ for what "means something" amounts to per field.

See the `Gridding`_ section of the AMReX documentation for details of how individual grids are created.

.. _`Gridding`: https://amrex-codes.github.io/amrex/docs_html/ManagingGridHierarchy_Chapter.html

Static Mesh Refinement
----------------------

For static refinement, we control the placement of grids by specifying
the low and high extents (in physical space) of each box in the lateral
directions.

The following example demonstrates how to tag regions for static refinement.
In this first example, all cells in the region :math:`[0.15,0.25,\texttt{prob_lo_z}] \times [0.35,0.45,\texttt{prob_hi_z}]`
and in the region :math:`[0.65,0.75,\texttt{prob_lo_z}]\times[0.85,0.95,\texttt{prob_hi_z}]` are tagged for
one level of refinement, where prob_lo_z and prob_hi_z are the vertical extents of the domain. They will be refined
by a factor of 2 in the x and y directions, and not refined (i.e. refinement ratio 1) in the z direction. If needed, the
refinement region will be expanded slightly to include the entirety of any partially-tagged coarse-level cell.

::

          amr.max_level = 1
          amr.ref_ratio = 2 2 1

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
          amr.ref_ratio = 3 3 1   9 9 1   #each triplet is refinement ratio in x,y,z for a single level

          remora.refinement_indicators = box1

          remora.box1.in_box_lo = .15 .25
          remora.box1.in_box_hi = .35 .45

And in this final example, the region :math:`[0.15,0.25,\texttt{prob_lo_z}]\times[0.35,0.45,\texttt{prob_hi_z}]`
will be refined by two levels of factor 3, but the larger region, :math:`[0.05,0.05,\texttt{prob_lo_z}]\times [0.75,0.75,\texttt{prob_hi_z}]``
will be refined by a single factor 3 refinement.

::

          amr.max_level = 2
          amr.ref_ratio = 3 3 1   9 9 1

          remora.refinement_indicators = box1 box2

          remora.box1.in_box_lo = .15 .25
          remora.box1.in_box_hi = .35 .45

          remora.box2.in_box_lo = .05 .05
          remora.box2.in_box_hi = .75 .75
          remora.box2.max_level = 1


We note that instead of specifying the physical extent enclosed, we can instead specify the indices of
the bounding box of the refined region in the index space of that fine level.
To do this we use
``in_box_lo_indices`` and ``in_box_hi_indices`` instead of ``in_box_lo`` and ``in_box_hi``.
If we want to refine the inner region (spanning half the width in each direction) by one level of
factor 2 refinement, and the domain has 32x64x8 cells at level 0 covering the domain, then we would set

::

          amr.max_level = 1
          amr.ref_ratio = 2 2 2

          remora.refinement_indicators = box1

          remora.box1.in_box_lo_indices = 16 32  4
          remora.box1.in_box_hi_indices = 47 95 11
          remora.box1.max_level = 1

There is also an option to specify the indices of the bounding box of the refined region in the index space of the coarser level, using
``in_box_lo_indices_crse`` and ``in_box_hi_indices_crse``.  This is useful when the user has a particular region in mind that they want to refine,
and they know the indices of that region on the coarser level but not on the finer level.  In this case, the code will automatically adjust the
indices to create a valid box at the finer level.

::

          amr.max_level = 1
          amr.ref_ratio = 2 2 2

          remora.refinement_indicators = box1

          remora.box1.in_box_lo_indices_crse = 16 32  4
          remora.box1.in_box_hi_indices_crse = 47 95 11
          remora.box1.max_level = 1


The lo_indices should be divisible by the refinement ratio, and the hi_indices should be one less than a number divisible by the refinement ratio.
There are no such requirements for the coarse level indices, since the code will adjust them as needed to create a valid box at the finer level.

Dynamic Mesh Refinement
-----------------------

Dynamically created tagging functions are based on runtime data specified in the inputs file.
These dynamically generated functions test on either state variables or derived variables
defined in REMORA_derive.cpp.

Available tests include

-  ``value_greater``: :math:`\text{field} \geq \text{threshold}`

-  ``value_less``: :math:`\text{field} \leq \text{threshold}`

-  ``adjacent_difference_greater``: :math:`\text{max}( | \text{difference between any nearest-neighbor cell} | ) \geq \text{threshold}`

The example below adds three user-named criteria:

- ``hi_temp``: cells with temperature greater than 10 on level 0, and greater than 20 on level 1 and higher. Triggers up to AMR level 3;
- ``tempdiff``: cells having a difference in temperature of 0.01 or more from that of any immediate neighbor. Triggers up to level 2, and only when the problem time is between 100 and 300 seconds;
- ``lo_vort``: cells with relative vorticity less than 0, and separately the region :math:`[0.25,0.25,\texttt{prob_lo_z}]\times[0.75,0.75,\texttt{prob_hi_z}]`.

Note that giving a field criterion an ``in_box_lo``/``in_box_hi`` does **not** restrict that criterion to
the box, as ``lo_vort`` might suggest: the box is refined, and the field test is applied over the
whole domain. Only a box-only indicator, with no ``field_name``, refines a region and nothing else.

Note that ``temp`` is the name of a state variable and ``vorticity`` is a derived variable.
Valid field options for refinement are any cell-centered tracer this run actually has, along with
``x_velocity``, ``y_velocity``, ``z_velocity``, ``vorticity``, ``mask``, and, in a build with particles,
``<particle>_count``. Which tracers exist depends on the input: ``temp`` and ``salt`` always, ``tracer``
and numbered ``tracer_1`` and up only when ``remora.nscalar`` asks for them (it defaults to 0), and
biology names such as ``NO3`` only with a biology model. Naming a field this run does not have aborts at
setup with the list of the ones it does. All but ``mask`` are restricted to water; see
`Masked Regions and Tagging`_ below for what that means field by field. Prefer ``value_greater`` for a particle count: its ghost cells are left at zero rather than
filled, so ``adjacent_difference_greater`` on one sees a step at every grid boundary and tags a set of
cells that depends on the domain decomposition.

::

          remora.refinement_indicators = hi_temp tempdiff lo_vort

          remora.hi_temp.max_level = 3
          remora.hi_temp.value_greater = 10. 20.
          remora.hi_temp.field_name = temp

          remora.tempdiff.max_level = 2
          remora.tempdiff.adjacent_difference_greater = 0.01
          remora.tempdiff.field_name = temp
          remora.tempdiff.start_time = 100
          remora.tempdiff.end_time = 300

          remora.lo_vort.max_level = 1
          remora.lo_vort.value_less = 0
          remora.lo_vort.field_name = vorticity
          remora.lo_vort.in_box_lo = .25 .25
          remora.lo_vort.in_box_hi = .75 .75

Masked Regions and Tagging
--------------------------

A field-based criterion -- ``value_greater``, ``value_less`` or ``adjacent_difference_greater`` -- on any
field but ``mask`` is evaluated only where the value it reads was computed from water alone, and
``adjacent_difference_greater`` differences two wet-cell values only. This mirrors how AMReX evaluates the
same criteria in the presence of an embedded boundary, where a covered cell is skipped and a difference
is taken only across a face the geometry leaves open.

Which cells that admits depends on where the field lives and how it is masked, so it is not simply "the
wet cells":

- a tracer (``temp``, ``salt``, ``tracer``, a biology tracer) is masked in place, so
  the test is just whether that cell is water;
- ``x_velocity`` and ``y_velocity`` are stored at a cell index but live on a face, and are masked by
  ``msku(i,j) = mskr(i-1,j) * mskr(i,j)`` and ``mskv(i,j) = mskr(i,j-1) * mskr(i,j)``. A water cell whose
  neighbor across that face is land therefore holds an exact zero that is a mask artifact rather than
  stationary water, so both cells sharing the face must be water;
- ``vorticity`` is a centered difference of the cell-centered velocities that is not itself masked, so
  its value depends on the whole 3x3 block of cells around it and all nine must be water. This is a
  workaround for the derived field being unmasked and can be narrowed once it is not;
- ``z_velocity`` is tested on its own cell, but only because nothing currently writes it: the vertical
  velocity the model solves for is held in a temporary, so the plotted and taggable ``z_velocity`` is
  identically zero. If it is ever connected, its dependence will be the five-point cross, not its own
  cell.

A land cell is never tagged by any of these. The one field exempt from all of it is ``mask``, below.

Nothing else is untagged. A region named explicitly with ``in_box_lo``/``in_box_hi`` (or the index-space
forms) is refined in full, land included, which is usually what is wanted when the region of interest
straddles a coast. Land may also end up refined because it is adjacent to a tagged region, or because
``amr.n_error_buf`` grew one.

To refine the coastline deliberately, use ``mask`` as the field name. It is the one field exempt from the
rule above -- a criterion keyed on the mask is asking where the coast is, so it is evaluated on every
cell. Since the mask is 0 or 1 exactly, ``adjacent_difference_greater = 0.5`` on it tags every cell whose
neighbor differs, which is the water cells along the coast and the land cells facing them, so the coast
is refined from both sides.

::

          remora.refinement_indicators = coast

          remora.coast.max_level = 1
          remora.coast.adjacent_difference_greater = 0.5
          remora.coast.field_name = mask

How the mask itself is carried across levels, and how it weights the two-way average, is described
in :ref:`Land/Sea Masking <sec:masking>`.

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

Example
-------

.. figure:: ./figures/scalar_whitebg_circle_hr_00010.png
    :width: 400

    Mesh refinement example for scalar advection. The black lines show the higher-resolution grids.

The ``Advection`` problem simulates the advection of a Gaussian-distributed passive scalar. The example above was generated with the inputs file found in ``Exec/Advection/inputs_ml``. The regions with scalar density greater than 0.5 are tagged for refinement after 200 seconds of evolution. In order to have the smaller level 1 refined grids shown above, it was run with the runtime parameter ``amr.max_grid_size=16 16 16``.

