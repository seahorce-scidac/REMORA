.. _Building:

Building
--------

REMORA can be built using either GNU Make or CMake. The following instructions apply to building on any system. We also provide instructions for building on :ref:`Perlmutter<Perlmutter (NERSC)>`.

Minimum Requirements
~~~~~~~~~~~~~~~~~~~~

REMORA requires a C++ compiler that supports the C++20 standard (required by AMReX) and a C compiler that supports the C99 standard.
Building with GPU support may be done with CUDA, HIP, or SYCL.
For CUDA, REMORA requires versions >= 12.2 (required by AMReX). For HIP and SYCL, only the latest compilers are supported.
Prerequisites for building with GNU Make include Python (>= 2.7, including 3) and standard tools available
in any Unix-like environments (e.g., Perl and sed). For building with CMake, the minimal requirement is version 3.25 (required by AMReX).

   .. note::
      **While REMORA is designed to work with SYCL, we do not make any guarantees that it will build and run on your Intel platform.**

Paradigm
~~~~~~~~

REMORA builds the executable in the ``Exec`` directory.
The problem directories within ``Exec`` include a number of problems which are used for testing.

.. _netcdf:

NetCDF (Optional)
~~~~~~~~~~~~~~~~~

REMORA uses `PnetCDF <https://parallel-netcdf.github.io/>`_ for optional NetCDF support. To build REMORA with PnetCDF, first install PnetCDF as per the instructions. Make a note of the directory where the library is installed, which we will call ``PNETCDF_DIR``. When compiling, add the PnetCDF ``pkgconfig`` directory to the environment variable ``PKG_CONFIG_PATH``, e.g.:

   .. code:: shell

             PKG_CONFIG_PATH=$PKG_CONFIG_PATH:$PNETCDF_DIR/lib/pkgconfig

At run-time, you may need to add PnetCDF to the link path, e.g.:

   .. code:: shell

             LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PNETCDF_DIR/lib

See sections below for compiler-specific instructions for how to enable NetCDF support.

   .. note::
      **When built with some MPI implementations, PnetCDF will crash when writing files larger than 2GB. Building
      with MPICH v4.2.2 will likely solve this issue.**

.. _netcdf4:

NetCDF-4 input files
^^^^^^^^^^^^^^^^^^^^

REMORA reads NetCDF-4 (HDF5) input files through PnetCDF's NetCDF-4 driver, which hands them to the netCDF-C library. The driver is off by default. To turn it on:

#. Install HDF5 with parallel support, and a netCDF-C built against it with parallel I/O (``nc-config --has-parallel4`` reports ``yes``). netCDF-C must not itself be built with PnetCDF support. Ubuntu, for example, provides both as ``libhdf5-mpi-dev`` and ``libnetcdf-mpi-dev``.
#. Configure PnetCDF with ``--with-netcdf4=<netCDF-C install directory>``. If netCDF-C's library has a name other than ``libnetcdf``, as Ubuntu's ``libnetcdf_mpi`` does, pass it with ``LIBS``, so that PnetCDF records it as a dependency, e.g.:

   .. code:: shell

             ./configure --prefix=$PNETCDF_DIR \
                         --with-netcdf4=/usr/lib/x86_64-linux-gnu/netcdf/mpi \
                         LIBS=-lnetcdf_mpi

   On Perlmutter, with ``cray-hdf5-parallel`` and ``cray-netcdf-hdf5parallel`` loaded (and ``cray-parallel-netcdf`` not loaded), the library is ``libnetcdf_parallel``:

   .. code:: shell

             ./configure --prefix=$PNETCDF_DIR --with-netcdf4=$NETCDF_DIR --disable-fortran \
                         CC=cc CXX=CC LIBS="-L$NETCDF_DIR/lib -lnetcdf_parallel"

   Without the right ``LIBS``, the build stops linking ``ncmpigen`` with undefined references to ``nc_open_par`` and the other ``nc_*`` functions.

   The configure summary should report ``NetCDF-4 support - enabled``, ``$PNETCDF_DIR/bin/pnetcdf-config --netcdf4`` should print ``enabled``, and ``ldd $PNETCDF_DIR/lib/libpnetcdf.so`` should list the netCDF-C and HDF5 libraries.
#. Check that ``pkg-config --libs pnetcdf`` works with ``$PNETCDF_DIR/lib/pkgconfig`` on ``PKG_CONFIG_PATH``. PnetCDF's ``pnetcdf.pc`` names ``netcdf`` as a requirement, so it fails with ``Package 'netcdf', required by 'pnetcdf', not found`` when netCDF-C has no ``netcdf.pc`` on the path, as on Perlmutter. Add netCDF-C's ``pkgconfig`` directory to ``PKG_CONFIG_PATH`` if it has one. Otherwise, since ``libpnetcdf.so`` already loads netCDF-C itself (the ``ldd`` check above), delete the requirement:

   .. code:: shell

             sed -i '/^[[:space:]]*Requires/d' $PNETCDF_DIR/lib/pkgconfig/pnetcdf.pc

#. Build REMORA against it as above, starting from an empty CMake cache if an earlier configure failed to find PnetCDF.

At run time, make sure the executable loads this PnetCDF library and not another one: ``LD_LIBRARY_PATH`` takes precedence over the library path recorded at build time, so a different PnetCDF on it is the one used. ``ldd remora_exec`` shows which one is loaded. A REMORA built against a PnetCDF without the driver, or run with one, stops with an error naming the NetCDF-4 file it could not open.

The ``NCFormats`` unit test (``ctest -R NCFormats``) writes and reads a small file in each format the linked PnetCDF supports, and reads a compressed NetCDF-4 file, so it checks that the driver works.


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

#. ``cd`` to  ``REMORA/Exec``

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
   |                 |                                  |                  |             |
   |                 | functionality enabled            |                  |             |
   +-----------------+----------------------------------+------------------+-------------+
   | USE_MOAB        | Whether to compile with MOAB     | TRUE / FALSE     | FALSE       |
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

   The name of the resulting executable (generated by the GNUmake system) encodes several of the build characteristics,
   including dimensionality of the problem, compiler name, and whether MPI and/or OpenMP were linked with the executable.
   Thus, several different build configurations may coexist simultaneously in a problem folder.
   For example, the executable may be named  ``REMORA3d.gnu.MPI.ex``, indicating that this is a 3d version of the code
   built with ``COMP=gnu``, and ``USE_MPI=TRUE``.


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

The options are described in the table below. Standard CMake options like ``CMAKE_INSTALL_PREFIX`` are also available.

   +----------------------------------+----------------------------------+------------------+-------------+
   | Option name                      | Description                      | Possible values  | Default     |
   |                                  |                                  |                  |             |
   |                                  |                                  |                  | value       |
   +==================================+==================================+==================+=============+
   | ``CMAKE_BUILD_TYPE``             | Build type                       | RelWithDebInfo / | Release     |
   |                                  |                                  |                  |             |
   |                                  |                                  | Release / Debug  |             |
   +----------------------------------+----------------------------------+------------------+-------------+
   | ``REMORA_ENABLE_MPI``            | Whether to enable MPI            | ON / OFF         | OFF         |
   +----------------------------------+----------------------------------+------------------+-------------+
   | ``REMORA_ENABLE_OPENMP``         | Whether to enable OpenMP         | ON / OFF         | OFF         |
   +----------------------------------+----------------------------------+------------------+-------------+
   | ``REMORA_ENABLE_CUDA``           | Whether to enable CUDA           | ON / OFF         | OFF         |
   +----------------------------------+----------------------------------+------------------+-------------+
   | ``REMORA_ENABLE_HIP``            | Whether to enable HIP            | ON / OFF         | OFF         |
   +----------------------------------+----------------------------------+------------------+-------------+
   | ``REMORA_ENABLE_SYCL``           | Whether to enable SYCL           | ON / OFF         | OFF         |
   +----------------------------------+----------------------------------+------------------+-------------+
   | ``REMORA_ENABLE_PNETCDF``        | Whether to compile with          | ON / OFF         | OFF         |
   |                                  |                                  |                  |             |
   |                                  | PNetCDF                          |                  |             |
   +----------------------------------+----------------------------------+------------------+-------------+
   | ``REMORA_ENABLE_PARTICLES``      | Whether to compile with          | ON / OFF         | OFF         |
   |                                  |                                  |                  |             |
   |                                  | particle functionality           |                  |             |
   |                                  |                                  |                  |             |
   |                                  | enabled                          |                  |             |
   +----------------------------------+----------------------------------+------------------+-------------+
   | ``REMORA_ENABLE_MOAB``           | Whether to build with MOAB       | ON / OFF         | OFF         |
   +----------------------------------+----------------------------------+------------------+-------------+
   | ``REMORA_ENABLE_TESTS``          | Whether to build tests           | ON / OFF         | OFF         |
   +----------------------------------+----------------------------------+------------------+-------------+
   | ``REMORA_ENABLE_FCOMPARE``       | Whether to build the AMReX       | ON / OFF         | OFF         |
   |                                  |                                  |                  |             |
   |                                  | utility fcompare                 |                  |             |
   +----------------------------------+----------------------------------+------------------+-------------+
   | ``REMORA_ENABLE_DOCUMENTATION``  | Whether to build                 | ON / OFF         | OFF         |
   |                                  |                                  |                  |             |
   |                                  | documentation                    |                  |             |
   +----------------------------------+----------------------------------+------------------+-------------+

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
