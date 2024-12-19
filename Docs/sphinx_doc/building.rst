.. _Building:

Building
--------

REMORA can be built using either GNU Make or CMake. The following instructions apply to building on any system. We also provide instructions for building on :ref:`Perlmutter<Perlmutter (NERSC)>`.

Minimum Requirements
~~~~~~~~~~~~~~~~~~~~

REMORA requires a C++ compiler that supports the C++17 standard and a C compiler that supports the C99 standard.
Building with GPU support may be done with CUDA, HIP, or SYCL.
For CUDA, REMORA requires versions >= 11.0. For HIP and SYCL, only the latest compilers are supported.
Prerequisites for building with GNU Make include Python (>= 2.7, including 3) and standard tools available
in any Unix-like environments (e.g., Perl and sed). For building with CMake, the minimal requirement is version 3.18.

   .. note::
      **While REMORA is designed to work with SYCL, we do not make any guarantees that it will build and run on your Intel platform.**

Paradigm
~~~~~~~~

REMORA uses the paradigm that different executables are built in different subdirectories within the ``Exec`` directory.  When
using gmake (see below), the user/developer should build in the directory of the selected problem.  When using
cmake (see below), separate executables are built for all of the problem directories listed in ``Exec/CMakeLists.txt``.
The problem directories within ``Exec`` include a number of problems, which are also used for ::.

NetCDF (Optional)
~~~~~~~~~~~~~~~~~

REMORA uses `PnetCDF <https://parallel-netcdf.github.io/>`_ for optional NetCDF support. To build REMORA with PnetCDF, first install PnetCDF as per the instructions. Make a note of the directory where the library is installed, which we will call ``PNETCDF_DIR``. When compiling, add the PnetCDF ``pkgconfig`` directory to the environment variable ``PKG_CONFIG_PATH``, e.g.:

   .. code:: shell

             PKG_CONFIG_PATH=$PKG_CONFIG_PATH:$PNETCDF_DIR/lib/pkgconfig

At run-time, you may need to add PnetCDF to the link path, e.g.:

   .. code:: shell

             LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PNETCDF_DIR/lib

See sections below for compiler-specific instructions for how to enable NetCDF support.

GNU Make
~~~~~~~~

The GNU Make system is best for use on large computing facility machines and production runs. With the GNU Make implementation, the build system will inspect the machine and use known compiler optimizations explicit to that machine if possible. These explicit settings are kept up-to-date by the AMReX project.

Using the GNU Make build system involves first setting environment variables for the directories of the dependencies of REMORA which is the repository of AMReX. AMReX is provided as a git submodule in REMORA and can be populated by using

   .. code:: shell

         git submodule init; git submodule update

in the REMORA repo, or before cloning by using
   .. code:: shell

         git clone --recursive <remora_repo>

Although submodules of these projects are provided, they can be placed externally as long as the ``<REPO_HOME>`` environment variables for each dependency is set correctly. An example of setting the ``<REPO_HOME>`` environment variables in the user's ``.bashrc`` is shown below:

::

   export REMORA_HOME=${HOME}/REMORA
   export AMREX_HOME=${REMORA_HOME}/Submodules/AMReX

The GNU Make system is set up to use the path to AMReX submodule by default, so it is not necessary to set
these paths explicitly, unless it is desired to do so. It is also possible to use an external version of
AMReX, downloaded by running

   .. code:: shell

             git clone https://github.com/amrex-codes/amrex.git

in which case the ``AMREX_HOME`` environment variable must point to the location where AMReX has been downloaded, which will take precedence over the default path to the submodule. If using bash shell,

::

   export AMREX_HOME=/path/to/external/amrex

or if using tcsh,

::

   setenv AMREX_HOME /path/to/external/amrex

#. ``cd`` to the desired build directory, e.g.  ``REMORA/Exec/Upwelling/``

#. Edit the ``GNUmakefile``; options include

   +-----------------+----------------------------------+------------------+-------------+
   | Option name     | Description                      | Possible values  | Default     |
   |                 |                                  |                  | value       |
   +=================+==================================+==================+=============+
   | COMP            | Compiler (gnu or intel)          | gnu / intel      | None        |
   +-----------------+----------------------------------+------------------+-------------+
   | USE_MPI         | Whether to enable MPI            | TRUE / FALSE     | FALSE       |
   +-----------------+----------------------------------+------------------+-------------+
   | USE_OMP         | Whether to enable OpenMP         | TRUE / FALSE     | FALSE       |
   +-----------------+----------------------------------+------------------+-------------+
   | USE_CUDA        | Whether to enable CUDA           | TRUE / FALSE     | FALSE       |
   +-----------------+----------------------------------+------------------+-------------+
   | USE_HIP         | Whether to enable HIP            | TRUE / FALSE     | FALSE       |
   +-----------------+----------------------------------+------------------+-------------+
   | USE_SYCL        | Whether to enable SYCL           | TRUE / FALSE     | FALSE       |
   +-----------------+----------------------------------+------------------+-------------+
   | DEBUG           | Whether to use DEBUG mode        | TRUE / FALSE     | FALSE       |
   +-----------------+----------------------------------+------------------+-------------+
   | USE_PNETCDF     | Whether to compile with PnetCDF  | TRUE / FALSE     | FALSE       |
   +-----------------+----------------------------------+------------------+-------------+
   | USE_PARTICLES   | Whether to compile with particle | TRUE / FALSE     | FALSE       |
   |                 | functionality enabled            |                  |             |
   +-----------------+----------------------------------+------------------+-------------+
   | PROFILE         | Include profiling info           | TRUE / FALSE     | FALSE       |
   +-----------------+----------------------------------+------------------+-------------+
   | TINY_PROFILE    | Include tiny profiling info      | TRUE / FALSE     | FALSE       |
   +-----------------+----------------------------------+------------------+-------------+
   | COMM_PROFILE    | Include comm profiling info      | TRUE / FALSE     | FALSE       |
   +-----------------+----------------------------------+------------------+-------------+
   | TRACE_PROFILE   | Include trace profiling info     | TRUE / FALSE     | FALSE       |
   +-----------------+----------------------------------+------------------+-------------+

   .. note::
      **Do not set both USE_OMP and USE_CUDA to true.**

   Information on using other compilers can be found in the AMReX documentation at
   https://amrex-codes.github.io/amrex/docs_html/BuildingAMReX.html .

#. Make the executable by typing

   .. code:: shell

      make

   The name of the resulting executable (generated by the GNUmake system) encodes several of the build characteristics, including dimensionality of the problem, compiler name, and whether MPI and/or OpenMP were linked with the executable.
   Thus, several different build configurations may coexist simultaneously in a problem folder.
   For example, the default build in ``REMORA/Exec/Upwelling`` will look
   like ``REMORA3d.gnu.MPI.ex``, indicating that this is a 3-d version of the code, made with
   ``COMP=gnu``, and ``USE_MPI=TRUE``.


Job info
~~~~~~~~

The build information can be accessed by typing

   .. code:: shell

      ./REMORA*ex --describe

in the directory where the executable has been built.


CMake
~~~~~

CMake is often preferred by developers of REMORA; CMake allows for building as well as easy testing and verification of REMORA through the use of CTest which is included in CMake. CTest functionality requires additional options, described in :ref:`Testing`.

Using CMake involves an additional configure step before using the ``make`` command. It is also expected that the user has cloned the REMORA repo with the ``--recursive`` option or performed ``git submodule init; git submodule update`` in the REMORA repo to populate its submodules.

To build with CMake, a user typically creates a ``build`` directory in the project directory and in that directory the ``cmake <options> ..`` command is used to configure the project before building it. REMORA provides an example build directory called ``Build`` with example scripts for performing the CMake configure. Once the CMake configure step is done, then the ``make`` command will build the executable.

An example CMake configure/build command to build REMORA without MPI. Replace the compilers with those installed on your system:

::

    cmake -DCMAKE_BUILD_TYPE:STRING=Release \
          -DREMORA_ENABLE_MPI:BOOL=OFF \
          -DCMAKE_CXX_COMPILER:STRING=g++ \
          .. && make

An example CMake configure/build command to build REMORA with MPI is listed below:

::

    cmake -DCMAKE_BUILD_TYPE:STRING=Release \
          -DREMORA_ENABLE_MPI:BOOL=ON \
          -DCMAKE_CXX_COMPILER:STRING=mpicxx \
          .. && make


An example CMake configure/build command to build REMORA with MPI, PnetCDF, and particles is listed below:

::

    cmake -DCMAKE_BUILD_TYPE:STRING=Release \
          -DREMORA_ENABLE_MPI:BOOL=ON \
          -DCMAKE_CXX_COMPILER:STRING=mpicxx \
          -DREMORA_ENABLE_PARTICLES:BOOL=ON \
          -DREMORA_ENABLE_PNETCDF:BOOL=ON \
          .. && make


Note that CMake is able to generate makefiles for the Ninja build system as well which will allow for faster building of the executable(s).

.. include:: perlmutter.rst
