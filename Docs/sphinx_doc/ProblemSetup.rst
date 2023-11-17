 .. role:: cpp(code)
    :language: c++

 .. _Problem:

Problem Setup
=============

Each problem setup with a different initial e.g. temperature profile and bathymetry has its own subdirectory within ``Exec``. To create a new problem, create a new subdirectory and copy the following files from another subdirectory:

* ``prob.cpp``
* ``prob.H``
* ``inputs``
* ``GNUmakefile`` (change the directory in ``ROMSX_PROBLEM_DIR`` to match)
* ``CMakeLists.txt`` (change the problem name to match)
* ``Make.package``
* ``amrvis.defaults`` (for visualization with AMRVis)

The file ``prob.cpp`` contains a number of functions that set the initial temeprature profile, as well as surface stress, mixing coefficients, and bathymetry. New problem-specific input parameters can be defined by adding a variable to the ``ProbParm`` class in ``prob.H``, and reading in the value in ``amrex_probinit`` in ``prob.cpp``. See the AMReX documentation on `ParmParse <https://amrex-codes.github.io/amrex/docs_html/Basics.html#parmparse>`_ for how to add parameters.
