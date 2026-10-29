.. _PerformanceTuning:

Performance tuning
==================

This page collects the build and run-time settings that decide how fast REMORA runs, and
says what each one does to the answers. Where a statement rests on a measurement, the
measurement is the Upwelling problem at 1024 x 2048 x 40 cells (84 million cells) on one
Perlmutter CPU node (AMD Milan, 128 cores) or one Perlmutter GPU node (four A100s), for
ten time steps. Scaling results for earlier releases are on the :ref:`Performance` page.

The short version:

* Build ``Release`` with the profiler on, and read the inclusive ``REMORA::Evolve()`` time.
* On CPUs, measure a few tile shapes (``remora.omp_tile_size``) for your problem. The tile
  shape changed the time per step by 2x in our tests, more than any other single setting.
* One MPI rank per core on a node is as fast as any hybrid mix on that node. Threads pay
  off at scale and in memory per node, not in single-node throughput.
* Leave floating-point contraction off (the default). It costs little and makes the
  answers independent of tiling, thread count and vectorization.
* Turn off the diagnostics you do not need: ``remora.sum_interval = -1``, ``remora.v = 0``,
  and plot and checkpoint intervals sized to what you will look at.

Build configuration
-------------------

**Build type.** ``CMAKE_BUILD_TYPE=Release`` (the CMake default) or ``DEBUG=FALSE`` with
GNU make. Debug builds are several times slower and change no answers.

**Profiler.** ``-DAMReX_TINY_PROFILE=ON`` (CMake) or ``TINY_PROFILE=TRUE`` (GNU make) adds
a small overhead and prints a table of timed regions at the end of the run. Two
tables appear: exclusive times (a region without the regions it calls) and inclusive
times. For comparisons between runs use the inclusive time of ``REMORA::Evolve()``; it
covers the time stepping and excludes initialization, the final plotfile and the profiler
itself. The exclusive table is the one to read when looking for the costliest kernel.

**Target CPU.** ``REMORA_ENABLE_NATIVE_ARCH=ON`` (CMake) or ``NATIVE_ARCH=TRUE`` (GNU make)
compiles for the CPU of the build machine (``-march=native``), which lets the stencil
kernels use that CPU's full vector width instead of the baseline x86-64 one. The binary
may then fail to start on an older CPU, which is why it is off by default; on a cluster
whose login and compute nodes differ, build on a compute node or use the compiler's
flag for the compute node's architecture instead. The gain depends on the problem and is
small when the step is limited by memory bandwidth (see `MPI ranks and OpenMP threads`_).
With contraction off (next item) it changes no answers: in our tests the ``-march=native``
build was bitwise identical to the default build.

.. _fp-contraction:

**Floating-point contraction.** Modern CPUs have a fused multiply-add instruction that
computes ``a*b + c`` with a single rounding, where the plain multiply and add round twice.
Compilers fuse ("contract") such pairs when the target has the instruction, which
``-march=native`` and the compiler wrappers of most clusters (Perlmutter's included)
enable. REMORA builds with contraction **off** by default: ``-ffp-contract=off`` for GCC,
Clang and the Intel LLVM compilers, ``-no-fma`` for classic Intel and ``-Mnofma`` for
NVHPC, applied to REMORA and to the AMReX sources of the same build.

The reason is reproducibility rather than accuracy. Which pairs a compiler fuses is a
heuristic, and it can fuse the same expression differently in the vectorized body of a
loop and in the scalar remainder that handles the last few cells of each row. Which cells
land in the remainder depends on where the tile and box edges fall, so with contraction
on, the answer depends on ``remora.omp_tile_size`` and on the domain decomposition. We
observed this in the vertical grid computation: the depths of the ``z`` levels differed in
their last bits between 8x8 tiles and 16-row strips, and the velocities differed at the
1e-8 level after five steps. With contraction off the arithmetic is exactly as written,
and the results are the same bitwise whatever the tiling, the thread count and the
vectorization (these are checked by the regression tests and the perf tests).

The cost is at most the fraction of the step spent in arithmetic that could have been
fused; on a bandwidth-bound step it is a few percent at most. Set
``REMORA_ENABLE_FP_CONTRACT=ON`` (CMake) or ``FP_CONTRACT=TRUE`` (GNU make) to allow
fusion again, with the tile-size dependence described above. Two further notes:

* The flag is a host-compiler flag and is applied to CPU builds only. GPU builds run no
  tiled host loops, so the issue does not arise there, and the device compilers keep their
  own defaults (``nvcc`` fuses by default). CPU and GPU results therefore differ at
  roundoff level, with or without this option.
* The regression tests are run without ``-march`` flags and so have never exercised
  contraction; a kernel that is tile-independent in ctest can still be tile-dependent in
  a contracting build. When changing kernels, also compare a ``-march=native`` build with
  contraction enabled at 8x8 tiles, 16-row strips and untiled.

MPI ranks and OpenMP threads
----------------------------

REMORA decomposes each level into boxes of at most ``amr.max_grid_size`` cells on a side
(``amr.max_grid_size_x`` etc. for different sizes per direction) and distributes the boxes
over the MPI ranks. Every rank should own at least one box; a rank with none idles. With
OpenMP, every ``MFIter`` loop of the time step is threaded over tiles (see `Tiling`_), so
the threads of a rank share its boxes. A rank whose box count times tile count is below
its thread count leaves threads idle; REMORA prints a warning at startup when a rank has
fewer boxes than threads and tiling is off.

On one Perlmutter CPU node with the 84-million-cell problem, the time per step at full
occupancy was within a few percent for 128 ranks x 1 thread, 64 x 2, 32 x 4 and 16 x 8
(each with 128 x 128 boxes and the tiling that suits the thread count). The node's memory
bandwidth sets the time once every core is busy, and ranks and threads draw on the same
bandwidth. Hybrid runs pay off elsewhere: fewer and larger messages and less halo volume
at scale, and less memory per node (ghost cells, buffers, duplicated metadata).

Thread placement matters. Bind each rank to a set of cores and spread its threads over
them, for example with Slurm ``srun -c <2 x threads> --cpu-bind=cores`` plus
``OMP_PROC_BIND=spread OMP_PLACES=cores`` (Perlmutter cores have two hardware threads; use
one thread per core). Open MPI's ``mpiexec`` binds a run of one or two ranks to a single
core by default, so pass ``--bind-to none`` or a ``--map-by`` option that reserves cores for
the threads; otherwise every thread of the rank shares one core.

Tiling
------

``remora.omp_tile_size`` (an alias of ``fabarray.mfiter_tile_size``) gives the maximum
tile extent in each direction; a box is split into tiles of that size, each ``MFIter``
loop visits the tiles, and with OpenMP the tiles of a rank are shared out among its
threads. The default on CPUs is ``1024000 16 1024000``, 16-row strips spanning the full
``x`` extent of a box; GPU builds run untiled. The third component may not be smaller than
the number of vertical levels.

The tile shape sets the working set of each kernel pass and so how much of it stays in
cache between passes. On the Milan node at one rank per core, column tiles of 16x16 and
8x8 cells ran 25 to 30 percent faster than the 16-row default, 8-row strips were in
between, and the untiled run took twice as long as the best. On the 4-core workstation
the code was developed on, the ordering of 16-row strips and 8x8 tiles was the reverse.
So the default is a reasonable starting point, not an answer: run a short timing of your
problem on your machine with

::

   remora.omp_tile_size = 8 8 1024000
   remora.omp_tile_size = 16 16 1024000
   remora.omp_tile_size = 1024000 16 1024000

and keep the fastest. When the box count is small relative to the thread count, smaller
tiles also mean more parallelism. The tile shape never changes the answers in the default
build (see :ref:`Floating-point contraction <fp-contraction>` for the exception).
``Tests/perf/README.md``
describes a script that runs such scans on Perlmutter and checks the plotfiles bitwise.

Grid layout and regridding
--------------------------

* ``amr.max_grid_size`` and ``amr.blocking_factor`` set the box sizes. One box per rank
  is the usual target on CPUs; on GPUs larger boxes (256 x 256 columns) amortize the
  kernel launches better and one box per GPU is best when the memory allows.
* Static nests do not need regridding. When the refined region never changes, set
  ``amr.regrid_int = -1`` so no time is spent tagging cells and rebuilding levels (see
  the notes on ``amr.regrid_int`` in :ref:`Inputs <sec:Inputs>`).
* The number of barotropic (fast) steps per baroclinic step, ``remora.ndtfast``, sets
  the cost of the 2D loop, which is the largest part of a single-level step. Use the
  smallest value that keeps the fast step stable.

Output and diagnostics
----------------------

* ``remora.sum_interval = -1`` turns off the integrated diagnostics (global sums of mass
  and tracers), which need a global reduction and a pass over the state each time they
  run.
* ``remora.v = 0`` suppresses the per-step reporting.
* ``remora.plot_int`` and ``remora.check_int`` control the plotfile and checkpoint
  frequencies. Writing a plotfile or checkpoint of a large domain costs the equivalent
  of several time steps, so size these to what you will look at and restrict
  ``remora.plot_vars_3d`` to the fields you need.

GPU
---

* Use one MPI rank per GPU, with ``amrex.use_gpu_aware_mpi = 1`` when the MPI library
  supports it, so halo exchanges move device memory directly.
* GPU builds run untiled; ``remora.omp_tile_size`` has no effect.
* The default box size gives a kernel launch per box per pass, so prefer few, large boxes
  per GPU. On the 84-million-cell problem with 256 x 256 boxes the four A100s of one
  Perlmutter node completed a step in a third of the time of the 128-core CPU node.
* Results differ from CPU results at roundoff level because the device compiler fuses
  multiply-adds (see :ref:`Floating-point contraction <fp-contraction>`).

Reproducibility
---------------

With the default build options the answers do not depend on the tile size, the number of
OpenMP threads or the compiler's vectorization choices; the regression suite has been
checked at one and two threads and at several tile sizes, and kernel changes are expected
to keep it that way. This relies on contraction being off and on two coding rules: a kernel
writes each cell of shared scratch from exactly one tile (or every tile that writes it
writes the same value), and sums whose order could vary, such as the coarse-fine flux
correction, are accumulated in a fixed order. Runs on different numbers of MPI ranks
agree bitwise for the same reasons as long as the box layout does not change the order
of such sums.
