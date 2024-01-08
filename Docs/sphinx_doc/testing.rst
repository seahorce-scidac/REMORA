.. _Testing:

Testing and Verification
------------------------

Testing and verification of REMORA can be performed using CTest, which is included in the CMake build system. If one builds REMORA with CMake, the testing suite, and the verification suite, can be enabled during the CMake configure step. A nightly test is reflected on the dashboard `here <https://my.cdash.org/index.php?project=REMORA>`_ .

An example ``cmake`` configure command performed in the ``Build`` directory in REMORA is shown below with options relevant to the testing suite:

::

  cmake -DCMAKE_INSTALL_PREFIX:PATH=./install \
        -DCMAKE_BUILD_TYPE:STRING=Release \
        -DREMORA_ENABLE_MPI:BOOL=ON \
        -DCMAKE_CXX_COMPILER:STRING=mpicxx \
        -DCMAKE_C_COMPILER:STRING=mpicc \
        -DCMAKE_Fortran_COMPILER:STRING=mpifort \
        -DREMORA_ENABLE_FCOMPARE:BOOL=ON \
        -DREMORA_ENABLE_TESTS:BOOL=ON \
        -DREMORA_USE_CPP:BOOL=ON \
        ..

While performing a ``cmake -LAH ..`` command will give descriptions of every option for the CMake project. Descriptions of particular options regarding the testing suite are listed below:

**REMORA_ENABLE_FCOMPARE** -- builds the ``fcompare`` utility from AMReX as well as the executable(s), to allow for testing differences between plot files

**REMORA_ENABLE_TESTS** -- enables the base level regression test suite that will check whether each test will run its executable to completion successfully


Building the Tests
~~~~~~~~~~~~~~~~~~

Once the user has performed the CMake configure step, the ``make`` command will build
every executable required for each test.
In this step, it is highly beneficial for the user to use the ``-j`` option for ``make``
to build source files in parallel.

Running the Tests
~~~~~~~~~~~~~~~~~

Once the test executables are built, CTest also creates working directories for each test within the ``Build`` directory
where plot files will be output, etc. This directory is analogous to the source location of the tests in ``Tests/test_files``.

To run the test suite, run ``ctest`` in the ``Build`` directory. CTest will run the tests and report their exit status.
Useful options for CTest are ``-VV`` which runs in a verbose mode where the output of each test can be seen. ``-R``
where a regex string can be used to run specific sets of tests. ``-j`` where CTest will bin pack and run tests in
parallel based on how many processes each test is specified to use and fit them into the amount of cores available
on the machine. ``-L`` where the subset of tests containing a particular label will be run.
(We note that using ``-L 'regression'`` will run all the tests that do not use SUNDIALS.)
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
regression test will refer to should be placed in ``Tests/REMORA-Gold-Files/<test_name>``. Next, edit the
``Exec/CMakeLists.txt`` and ``Tests/CTestList.cmake`` files, add the problem and the corresponding tests
to the list. Note that there are different categories of tests and if your test falls outside of these
categories, a new function to add the test will need to be created. After these steps, your test will be
automatically added to the test suite database when doing the CMake configure with the testing suite enabled.
