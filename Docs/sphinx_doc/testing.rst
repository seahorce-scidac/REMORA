.. _Testing:

Testing and Verification
------------------------

Testing and verification of REMORA can be performed using CTest, which is included in the CMake build system. If one builds REMORA with CMake, the testing suite, and the verification suite, can be enabled during the CMake configure step.

An example ``cmake`` configure/build command performed in the ``Build`` directory in REMORA is shown below with options relevant to the testing suite:

::

  cmake -DCMAKE_INSTALL_PREFIX:PATH=./install \
        -DCMAKE_BUILD_TYPE:STRING=Release \
        -DREMORA_ENABLE_MPI:BOOL=ON \
        -DCMAKE_CXX_COMPILER:STRING=mpicxx \
        -DREMORA_ENABLE_FCOMPARE:BOOL=ON \
        -DREMORA_ENABLE_TESTS:BOOL=ON \
        -DREMORA_USE_CPP:BOOL=ON \
        .. && make

To run CTest on GPU, add one of the options: ``REMORA_ENABLE_CUDA``, ``REMORA_ENABLE_HIP``, or ``REMORA_ENABLE_SYCL``, depending on the type of GPU on your system. For example:

::

  cmake -DCMAKE_INSTALL_PREFIX:PATH=./install \
        -DCMAKE_BUILD_TYPE:STRING=Release \
        -DREMORA_ENABLE_MPI:BOOL=ON \
        -DREMORA_ENABLE_CUDA:BOOL=ON \
        -DCMAKE_CXX_COMPILER:STRING=mpicxx \
        -DREMORA_ENABLE_FCOMPARE:BOOL=ON \
        -DREMORA_ENABLE_TESTS:BOOL=ON \
        -DREMORA_USE_CPP:BOOL=ON \
        .. && make

While performing a ``cmake -LAH ..`` command will give descriptions of every option for the CMake project. Descriptions of particular options regarding the testing suite are listed below:

**REMORA_ENABLE_FCOMPARE** -- builds the ``fcompare`` utility from AMReX as well as the executable(s), to allow for testing differences between plot files

**REMORA_ENABLE_TESTS** -- enables the base level regression test suite that will check whether each test will run its executable to completion successfully, and the unit tests in ``Tests/Unit``


Building the Tests
~~~~~~~~~~~~~~~~~~

Once the user has performed the CMake configure step, the ``make`` command will build
every executable required for each test.
In this step, it is highly beneficial for the user to use the ``-j`` option for ``make``
to build source files in parallel. For example:

   .. code:: shell

        make -j8

Running the Tests
~~~~~~~~~~~~~~~~~

Once the test executables are built, CTest also creates working directories for each test within the ``Build`` directory
where plot files will be output, etc. This directory is analogous to the source location of the tests in ``Tests/test_files``.

To run the test suite, run ``ctest`` in the ``Build`` directory. CTest will run the tests and report their exit status.
Useful options for CTest are ``-VV`` which runs in a verbose mode where the output of each test can be seen. ``-R``
where a regex string can be used to run specific sets of tests. ``-j`` where CTest will bin pack and run tests in
parallel based on how many processes each test is specified to use and fit them into the amount of cores available
on the machine.
Output for the last set of tests run is available in the ``Build`` directory in ``Testing/Temporary/LastTest.log``.

Adding Tests
~~~~~~~~~~~~

Developers are encouraged to add tests to REMORA and in this section we describe how the tests are organized in the
CTest framework. The locations (relative to the REMORA code base) of the tests are in ``Tests``. To add a test, first
create a problem directory with a name in ``Exec/<prob_name>``. This problem directory is meant for a production
run where the simulation is run until convergence or a solution is developed. This problem setup could comprise
of a more complex physics than the corresponding tests for regression at ``Tests/test_files/<test_name>``. Prepare
toned down versions of the input file(s) for each combination of physics that a regression test is desired.
For example, ``Upwelling`` problem with input file ``Exec/Upwelling/inputs`` solves the double gyre problem. The corresponding regression test is driven by the input files
``Tests/test_files/Upwelling/Upwelling.i``.

Any file in the test directory will be copied during CMake configure to the test's working directory.
The input files meant for regression test run only until a few time steps. The reference solution that the
regression test will refer to should be placed in ``Tests/REMORA_Gold_Files/<test_name>``. Next, edit the
``Exec/CMakeLists.txt`` and ``Tests/CTestList.cmake`` files, add the problem and the corresponding tests
to the list. Note that there are different categories of tests and if your test falls outside of these
categories, a new function to add the test will need to be created. After these steps, your test will be
automatically added to the test suite database when doing the CMake configure with the testing suite enabled.

Because the copy happens at CMake *configure* time via ``file(GLOB)``, adding or changing a file in a test
directory does not reach the working directory on an incremental ``cmake --build``. Re-run ``cmake`` after
adding a test or editing its input.

Kinds of test
~~~~~~~~~~~~~

The functions in ``Tests/CTestList.cmake`` cover several kinds of assertion. A snapshot against a gold file
is the default, but it only says "this run still does what it did when the baseline was blessed" -- it cannot
say the baseline was right, and it cannot notice a feature that has silently stopped doing anything.

``add_test_r``, ``add_test_r_hitol``
    Snapshot against ``Tests/REMORA_Gold_Files/<test_name>``, at 1e-11 and 1e-5 respectively.

``add_test_r_gold``
    Snapshot against *another* test's gold file. Use it when a run is expected to reproduce an existing
    baseline exactly -- for instance a high-resolution-bathymetry lane over a problem whose bathymetry is
    constant, where the average-down has to be a no-op. Costs no new gold data.

``add_test_r_differ``
    Must agree with its own gold **and** disagree with another one. This is the shape to use for a feature
    whose purpose is to change the answer: a lane that stops taking effect (a misspelled parameter, a dropped
    branch) still matches its own snapshot, and only the disagreement clause catches it. Both clauses are
    needed, since ``fcompare`` aborts outright on a level-count mismatch and a bare disagreement test would
    then pass for the wrong reason.

``add_test_equiv``
    Two inputs in one test directory, describing the same configuration by different routes, must agree.
    Neither run is a baseline, so a translation layer is covered without a gold file blessed by the code
    under test, and neither route can drift alone. The second input needs a different plotfile prefix,
    since both runs share a working directory. Agreement is blind to a value that both routes read wrongly
    in the same way -- two runs that each ignore an input agree perfectly -- so trailing
    ``<tol> <var> <min> <max> ...`` arguments add a ``check_extrema.sh`` assertion on the first run's
    plotfile, pinning the magnitude of whatever must not go degenerate.

``add_test_extrema``
    Assert a variable's min and max against values known from outside REMORA -- a closed-form reference, or
    a constant an initial condition must reproduce. Takes any number of ``<var> <min> <max>`` triples after
    the tolerance, all checked against one model run. Unlike a gold file this says what the number should
    *be*, so it also catches a baseline that was wrong when it was blessed, and it needs no baseline bytes.

``add_test_abort``
    Assert that a misconfigured input fails, and fails for the stated reason. The run must exit nonzero and
    the log must contain the given plain-text substring, so a successful run, a different abort, and a
    segfault all fail the test. Prefer this to ``WILL_FAIL``, which any nonzero exit satisfies.

``add_test_r_selfcompare``
    Run one input twice under two sets of command-line options and require the plotfiles to agree to
    1e-11. For an invariant between two code paths, such as a driver that must reproduce another at a
    degenerate setting, or output frequency not changing the answer. Needs no gold file.

``add_test_restart``
    Run an input to the end with a checkpoint partway, restart a second run from the named checkpoint,
    and require the two final plotfiles to agree: a restart must continue the run it was saved from.

``add_test_0``
    Require the final plotfile to equal ``plt00000`` to 1e-14. For a case whose exact solution is rest.

``add_test_conservation``
    Assert how far an integrated quantity drifts over a run, from the integrated-quantity log. The second
    argument names the input file, so several tests can share one case under different options, and the
    trailing ``<column> <bound> <below|above>`` triples are checked by ``check_conservation.sh``. An
    ``above`` bound is a control: it asserts the case drifts without the correction under test, so its
    ``below`` partner cannot pass because the case has stopped exercising anything.

``add_test_log``
    The run must succeed **and** its log must contain the given substring. For a code path whose answers
    are not worth a gold file but which must keep reaching a behavior it reports, or for a diagnostic that
    aborts on failure and prints a line on success.

``add_test_nlevels``, ``add_test_extent``, ``add_test_cells``
    Assert the refined grid: the deepest level that exists, the physical extent a level spans in one
    direction, or the number of cells a level covers, summed over boxes. A tagging bug still writes a
    self-consistent plotfile, so grid structure is the only place it shows. The first two also require
    the run to stay at rest, as ``add_test_0`` does; ``add_test_cells`` does not, for a case that evolves.
    Pin ``amr.n_error_buf``, ``amr.blocking_factor`` and ``amr.grid_eff`` in the input: each moves the grids
    without any tagging having changed, and ``main.cpp`` sets the first two rather than taking AMReX's
    defaults.

Cheap lanes are worth adding: ``remora.max_step = 0`` with ``remora.plot_int = 1`` is a legal
initialize-and-dump run of about a second, which is enough for ``add_test_extrema`` to pin an initial
condition with no time stepper in the way.

Unit tests
~~~~~~~~~~

Unit tests that exercise individual functions live in ``Tests/Unit``, built as ``remora_unit_tests`` and
labelled ``unit`` so ``ctest -L unit`` runs them alone and ``ctest -LE unit`` skips them.

There is no framework: a unit test is a ``main()`` that checks answers, prints a line per failure, and
returns nonzero. To add one, add a ``.cpp``
beside the existing test and list it in ``Tests/Unit/CMakeLists.txt``. The target links AMReX, so a test
sees the same precision as the build around it; a case that only holds in double has to be gated on that
rather than assumed.

The unit test ``REMORA_DateClock_test.cpp`` tests that the date and calendar functions in
``REMORA_DateClock.H`` match ROMS's ``dateclock.F``.

Regenerating gold files
~~~~~~~~~~~~~~~~~~~~~~~

A gold file is only meaningful if what produced it is known, so when regenerating one:

#. Work from a clean tree. ``git status --porcelain`` must be empty -- a plotfile's ``job_info`` records the
   REMORA git hash, and a ``-dirty`` hash makes the baseline unreproducible.
#. Delete the build directory first, or ``AMReX_buildInfo.cpp`` may be stale and record a build date and
   hash from an older configure.
#. Record what is not in ``job_info``: the number of MPI ranks used, and the SHA-256 of every input file
   the run read, including NetCDF data. Gold files are shared between the serial and 2-rank paths, so a
   baseline is only valid if the run is rank-invariant to the accepted tolerance.
#. Treat a change to any input datum as invalidating the baseline: the same commit should update the data,
   the gold file, and the recorded checksums together.
