 .. role:: cpp(code)
    :language: c++

.. _ProblemSetup:

Problem Setup
=============

Each problem setup with a different initial e.g. temperature profile and bathymetry has its own subdirectory within ``Exec``.
To create a new problem, create a new subdirectory and copy the following files from another subdirectory:

* ``inputs``
* ``amrvis.defaults`` (for visualization with AMRVis)

You may also copy the ``Exec/BlankProblem`` directory.

The files in ``Source/Prob`` contain a number of functions that set initial conditions (temperature, salinity, scalar, velocities),
as well as other variables like bathymetry, Coriolis forcing, and wind speed).  Which of these routines is called is controlled
by ``Exec/REMORA_prob.cpp``; if you introduce a new type of problem you will need to edit ``Exec/REMORA_prob.cpp`` and potentially add
new files in ``Source/Prob``.

New problem-specific input parameters can be defined by reading in additional values with the ParmParse prefix ``remora.prob``.
See the AMReX documentation on `ParmParse <https://amrex-codes.github.io/amrex/docs_html/Basics.html#parmparse>`_ for how to add parameters.

REMORA will call ``FillBoundary`` on the relevant variables after the functions called from ``REMORA_Prob.cpp``.
This will fill the values in the `ghost cells <https://amrex-codes.github.io/amrex/docs_html/Basics.html#ghost-cells>`_ at interior grid-grid and periodic boundaries.
However, it will not do so at interior grid-grid boundaries that fall on non-periodic domain boundaries.
This is primarily a concern for variables such as bathymetry, which are not specified by boundary conditions, but still have well-defined values at boundaries.

By default, these functions will be evaluated at level 0 and interpolated to higher levels as needed for mesh refinement. If ``remora.hires_grid_level > 0``, the analytic function for bathymetry is evaluated on level ``remora.hires_grid_level``. Bathymetry on higher levels is interpolated as usual. Bathymetry on lower levels is calculated by averaging down.
Likewise, if ``remora.hires_init_level > 0``, the initial-condition, sea surface height, and biology functions are evaluated on that level and averaged down.

That level does not exist yet when this happens, so the initial-condition and sea surface height functions must read coordinates (``z_r``, ``z_w``, ``x_r``, ``y_r``)
from the ``ProbCoords`` argument, which is laid out like the ``MultiFab`` being filled, rather than from the level arrays on the ``REMORA`` object.
