 .. role:: cpp(code)
    :language: c++

.. _GettingStarted:

Getting Started
===============

Quickstart
----------

This Quickstart guide will guide the user through downloading the code and building/running an Advection test case with GNUMake without MPI or OpenMP.
For building with GNUMake, REMORA requires a C++ compiler that supports the C++17 standard and a C compiler that supports the C99 standard.
The code is available on Github and can be accessed with ``git``.

   .. code:: shell

      git clone --recursive https://github.com/seahorce-scidac/REMORA.git

Now we enter the directory for the problem we're interested in, in this case Advection:

   .. code:: shell

      cd REMORA/Exec/Advection

And build, which will produce the executable ``REMORA.3d.gnu.TEST.ex``.

   .. code:: shell

      make -j USE_MPI=FALSE

The executable can be run with the provided inputs file:

   .. code:: shell

      ./REMORA.3d.gnu.TEST.ex inputs

This will produce an AMReX plotfile at the 10th time step called ``plt00010`` which can be :ref:`visualized<Visualization>`.

A similar process can be used to build other cases within ``Exec``, except for ``IdealMiniGrid``, which requires MPI and :ref:`PnetCDF<netcdf>`.

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
