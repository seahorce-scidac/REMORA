 .. role:: cpp(code)
    :language: c++

.. _GettingStarted:

Getting Started
===============

Quickstart
----------

This Quickstart guide will guide the user through downloading the code and building/running an Advection test case with GNUMake with MPI.
For building with GNUMake, REMORA requires a C++ compiler that supports the C++20 standard and a C compiler that supports the C99 standard.
The code is available on Github and can be accessed with ``git``.

   .. code:: shell

      git clone --recursive https://github.com/seahorce-scidac/REMORA.git

Now enter the Exec directory

   .. code:: shell

      cd REMORA/Exec

And build,

   .. code:: shell

      make -j USE_MPI=TRUE

which will produce the executable ``REMORA.3d.gnu.TEST.MPI.ex``.

Now to run the Advection problem, for example,

   .. code:: shell

      cd Advection

Then, for a single-rank run,

   .. code:: shell

      ../REMORA.3d.gnu.TEST.MPI.ex inputs

This will produce an AMReX plotfile at the 10th time step called ``plt00010`` which can be :ref:`visualized<Visualization>`.

A similar process can be used to build other cases within ``Exec``, except for ``IdealMiniGrid``, which requires :ref:`PnetCDF<netcdf>`.

Note, to build with PnetCDF support, set ``USE_PNETCDF=TRUE`` when invoking ``make``. This can be installed with spack following the instructions in the :ref:`Building<building>` section.


Downloading the code
--------------------

First, make sure that git is installed on your machine.

Then download the REMORA repository by typing:

   .. code:: shell

             git clone https://github.com/seahorce-scidac/REMORA.git

Or, to automatically include the AMReX submodule when downloading REMORA,
type:

   .. code:: shell

             git clone --recursive https://github.com/seahorce-scidac/REMORA.git

.. include:: submodule.rst

.. include:: building.rst

.. include:: InputFiles.rst


Building parallel-netCDF support with spack
-------------------------------------------

Clone Spack from GitHub into your home directory:

.. code:: shell

   cd ~
   git clone -c feature.manyFiles=true --depth=2 https://github.com/spack/spack.git

Set up the environment. The following command is for the Fish shell:

.. code:: shell

   source ~/spack/share/spack/setup-env.fish

Optionally, create a Spack environment:

.. code:: shell

   spack env create remora

Find and register the available compilers:

.. code:: shell

   spack compiler find

In the macOS setup used for these instructions, the compiler configuration was
stored in ``~/.spack/darwin/compilers.yaml``.

If bootstrapping is needed before installing packages, run:

.. code:: shell

   spack bootstrap now

Install parallel-netCDF and Open MPI:

.. code:: shell

   spack install parallel-netcdf
   spack install openmpi

.. note::

   You can specify a compiler explicitly, for example:

   .. code:: shell

      spack install openmpi%gcc
      spack install parallel-netcdf%gcc

Load the installed packages:

.. code:: shell

   spack load parallel-netcdf
   spack load openmpi
