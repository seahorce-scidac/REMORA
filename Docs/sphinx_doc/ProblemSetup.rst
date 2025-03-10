 .. role:: cpp(code)
    :language: c++

.. _ProblemSetup:

Problem Setup
=============

Each problem setup with a different initial e.g. temperature profile and bathymetry has its own subdirectory within ``Exec``. To create a new problem, create a new subdirectory and copy the following files from another subdirectory:

* ``prob.cpp``
* ``prob.H``
* ``inputs``
* ``GNUmakefile`` (change the directory in ``REMORA_PROBLEM_DIR`` to match)
* ``CMakeLists.txt`` (change the problem name to match)
* ``Make.package``
* ``amrvis.defaults`` (for visualization with AMRVis)

You may also copy the ``Exec/BlankProblem`` directory. The problem names will still need to be changed in ``GNUMakefile`` and ``CMakeLists.txt``.

The file ``prob.cpp`` contains a number of functions that set initial conditions (temperature, salinity, scalar, velocities), as well as other variables like bathymetry, Coriolis forcing, and wind speed). Examples of these functions can be found in other problem directories within ``Exec``. If they are not being used (i.e. the quantity is being specified by file rather than analytically), they do not need to appear in ``prob.cpp`` and ``prob.H``. The MultiFabs of the variables to be modified are passed in. Additional state variables can be accessed via the instance of the ``REMORA`` class that is passed in. For example, to get read-only access to an ``Array4`` of ``z_w`` within a MFIter, an example of which can be found in ``Exec/Upwelling/prob.cpp``:

.. code::

    Array4<const Real> const& z_w = remora.vec_z_w[lev]->const_array(mfi);

The current time can also be accessed using the method given below. An example is also found in ``Exec/Upwelling/prob.cpp``:

.. code::

    Real time = remora.get_t_old(lev);

New problem-specific input parameters can be defined by adding a variable to the ``ProbParm`` class in ``prob.H``, and reading in the value in ``amrex_probinit`` in ``prob.cpp``. See the AMReX documentation on `ParmParse <https://amrex-codes.github.io/amrex/docs_html/Basics.html#parmparse>`_ for how to add parameters.

REMORA will call ``FillBoundary`` on the relevant variables after the functions in ``prob.cpp`` are called. This will fill the values in the `ghost cells <https://amrex-codes.github.io/amrex/docs_html/Basics.html#ghost-cells>`_ at interior grid-grid and periodic boundaries. However, it will not do so at interior grid-grid boundaries that fall on non-periodic domain boundaries. This is primarily a concern for variables such as bathymetry, which are not specified by boundary conditions, but still have well-defined values at boundaries.
