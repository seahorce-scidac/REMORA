
# Have CMake discover the number of cores on the node
include(ProcessorCount)
ProcessorCount(PROCESSES)

set(FCOMPARE_GOLD_FILES_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/REMORA_Gold_Files)

#=============================================================================
# Functions for adding tests / Categories of tests
#=============================================================================
function(resolve_test_exe TEST_DIR TEST_EXE OUT_VAR)
    if(WIN32)
        # Multi-config generators place binaries in a config subdir.
        set(${OUT_VAR} "${CMAKE_BINARY_DIR}/Exec/${TEST_DIR}/*/${TEST_EXE}.exe" PARENT_SCOPE)
    else()
        set(_exe_in_subdir "${CMAKE_BINARY_DIR}/Exec/${TEST_DIR}/${TEST_EXE}${CMAKE_EXECUTABLE_SUFFIX}")
        set(_exe_in_root  "${CMAKE_BINARY_DIR}/Exec/${TEST_EXE}${CMAKE_EXECUTABLE_SUFFIX}")
        if(EXISTS "${_exe_in_subdir}")
            set(${OUT_VAR} "${_exe_in_subdir}" PARENT_SCOPE)
        elseif(EXISTS "${_exe_in_root}")
            set(${OUT_VAR} "${_exe_in_root}" PARENT_SCOPE)
        else()
            # Keep the historical path so the error message is still informative.
            set(${OUT_VAR} "${_exe_in_subdir}" PARENT_SCOPE)
        endif()
    endif()
endfunction()

macro(setup_test)
    set(CURRENT_TEST_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${TEST_NAME})
    set(CURRENT_TEST_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/test_files/${TEST_NAME})
    set(PLOT_GOLD ${FCOMPARE_GOLD_FILES_DIRECTORY}/${TEST_NAME})

    file(MAKE_DIRECTORY ${CURRENT_TEST_BINARY_DIR})
    file(GLOB TEST_FILES "${CURRENT_TEST_SOURCE_DIR}/*")
    file(COPY ${TEST_FILES} DESTINATION "${CURRENT_TEST_BINARY_DIR}/")

    if(REMORA_ENABLE_MPI)
        set(NP 2)
        set(MPI_COMMANDS "${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${NP} ${MPIEXEC_PREFLAGS}")
    else()
        set(NP 1)
        unset(MPI_COMMANDS)
    endif()

    # Set some default runtime options for all tests in this category
    # set(RUNTIME_OPTIONS "time.max_step=10 amr.plot_file=plt time.plot_interval=10 amrex.throw_exception=1 amrex.signal_handling=0")
    # set(RUNTIME_OPTIONS "max_step=10 amr.plot_file=plt amr.checkpoint_files_output=0 amr.plot_files_output=1 amrex.signal_handling=0")

endmacro(setup_test)

# Standard regression test
function(add_test_r TEST_NAME TEST_EXE PLTFILE)

    setup_test()

    resolve_test_exe("${TEST_DIR}" "${TEST_EXE}" TEST_EXE)

    set(FCOMPARE_TOLERANCE "-r 1e-11 --abs_tol 1.0e-11")
    set(FCOMPARE_FLAGS "-a ${FCOMPARE_TOLERANCE}")
    set(test_command sh -c "${MPI_COMMANDS} ${TEST_EXE} ${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.i ${RUNTIME_OPTIONS} > ${TEST_NAME}.log && ${FCOMPARE_EXE} ${FCOMPARE_FLAGS} ${PLOT_GOLD} ${CURRENT_TEST_BINARY_DIR}/${PLTFILE}")

    add_test(${TEST_NAME} ${test_command})
    set_tests_properties(${TEST_NAME}
        PROPERTIES
        TIMEOUT 5400
        PROCESSORS ${NP}
        WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
        LABELS "regression"
        ATTACHED_FILES_ON_FAIL "${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.log"
    )
endfunction(add_test_r)

function(add_test_r_hitol TEST_NAME TEST_EXE PLTFILE)
    setup_test()

    resolve_test_exe("${TEST_DIR}" "${TEST_EXE}" TEST_EXE)

    set(FCOMPARE_TOLERANCE "-r 1e-5 --abs_tol 1.0e-5")
    set(FCOMPARE_FLAGS "-a ${FCOMPARE_TOLERANCE}")
    set(test_command sh -c "${MPI_COMMANDS} ${TEST_EXE} ${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.i ${RUNTIME_OPTIONS} > ${TEST_NAME}.log && ${FCOMPARE_EXE} ${FCOMPARE_FLAGS} ${PLOT_GOLD} ${CURRENT_TEST_BINARY_DIR}/${PLTFILE}")

    add_test(${TEST_NAME} ${test_command})
    set_tests_properties(${TEST_NAME}
        PROPERTIES
        TIMEOUT 5400
        PROCESSORS ${NP}
        WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
        LABELS "regression"
        ATTACHED_FILES_ON_FAIL "${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.log"
    )
endfunction(add_test_r_hitol)

# Regression test comparing against ANOTHER test's gold file. Use when a run is expected to
# reproduce an existing baseline exactly -- e.g. a high-resolution-bathymetry lane over a
# problem whose bathymetry is constant, where the average-down must be a no-op. Costs no new
# gold data, which is the point.
function(add_test_r_gold TEST_NAME TEST_EXE PLTFILE GOLD_NAME)

    setup_test()

    resolve_test_exe("${TEST_DIR}" "${TEST_EXE}" TEST_EXE)

    set(FCOMPARE_TOLERANCE "-r 1e-11 --abs_tol 1.0e-11")
    set(FCOMPARE_FLAGS "-a ${FCOMPARE_TOLERANCE}")
    set(test_command sh -c "${MPI_COMMANDS} ${TEST_EXE} ${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.i > ${TEST_NAME}.log && ${FCOMPARE_EXE} ${FCOMPARE_FLAGS} ${FCOMPARE_GOLD_FILES_DIRECTORY}/${GOLD_NAME} ${CURRENT_TEST_BINARY_DIR}/${PLTFILE}")

    add_test(${TEST_NAME} ${test_command})
    set_tests_properties(${TEST_NAME}
        PROPERTIES
        TIMEOUT 5400
        PROCESSORS ${NP}
        WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
        LABELS "regression"
        ATTACHED_FILES_ON_FAIL "${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.log"
    )
endfunction(add_test_r_gold)

# Regression test that must agree with its OWN gold and DISAGREE with another one. Use for a
# feature whose whole point is to change the answer: a lane that silently stops doing anything
# (a misspelled parameter, a dropped branch) still matches its own gold, and only the
# disagreement clause catches it. Both clauses are needed -- fcompare aborts outright on a
# level-count mismatch, so a bare disagreement test would pass for the wrong reason.
function(add_test_r_differ TEST_NAME TEST_EXE PLTFILE OTHER_GOLD)

    setup_test()

    resolve_test_exe("${TEST_DIR}" "${TEST_EXE}" TEST_EXE)

    set(FCOMPARE_TOLERANCE "-r 1e-11 --abs_tol 1.0e-11")
    set(FCOMPARE_FLAGS "-a ${FCOMPARE_TOLERANCE}")
    set(test_command sh -c "${MPI_COMMANDS} ${TEST_EXE} ${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.i > ${TEST_NAME}.log && ${FCOMPARE_EXE} ${FCOMPARE_FLAGS} ${PLOT_GOLD} ${CURRENT_TEST_BINARY_DIR}/${PLTFILE} && ! ${FCOMPARE_EXE} ${FCOMPARE_FLAGS} ${FCOMPARE_GOLD_FILES_DIRECTORY}/${OTHER_GOLD} ${CURRENT_TEST_BINARY_DIR}/${PLTFILE}")

    add_test(${TEST_NAME} ${test_command})
    set_tests_properties(${TEST_NAME}
        PROPERTIES
        TIMEOUT 5400
        PROCESSORS ${NP}
        WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
        LABELS "regression"
        ATTACHED_FILES_ON_FAIL "${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.log"
    )
endfunction(add_test_r_differ)

# Test that a misconfigured input aborts, and aborts for the stated reason. The run must exit
# nonzero AND the log must carry the message, so a successful run, a different abort, and a
# segfault all fail. Deliberately not WILL_FAIL (which any nonzero exit satisfies) and not
# PASS_REGULAR_EXPRESSION (the model's stdout is redirected into the log, away from CTest).
# ABORT_SUBSTRING must be plain text. Note there is no ";" in the command: CMake would treat it
# as a list separator and hand sh -c two arguments, silently dropping everything after it.
function(add_test_abort TEST_NAME TEST_EXE ABORT_SUBSTRING)

    setup_test()

    resolve_test_exe("${TEST_DIR}" "${TEST_EXE}" TEST_EXE)

    set(test_command sh -c "! ${MPI_COMMANDS} ${TEST_EXE} ${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.i > ${TEST_NAME}.log 2>&1 && grep -q -- \"${ABORT_SUBSTRING}\" ${TEST_NAME}.log")

    add_test(${TEST_NAME} ${test_command})
    set_tests_properties(${TEST_NAME}
        PROPERTIES
        TIMEOUT 600
        PROCESSORS ${NP}
        WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
        LABELS "regression"
        ATTACHED_FILES_ON_FAIL "${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.log"
    )
endfunction(add_test_abort)

# Assert variables' min and max in a plotfile against values known from outside REMORA -- a
# closed-form reference, or constants that an initial condition must reproduce. Unlike a gold
# file this says what the numbers should BE, so it also catches a baseline that was wrong when
# it was blessed. Pass any number of <var> <min> <max> triples after TOL; they are all checked
# against a single model run. See Tests/check_extrema.sh for the comparison itself.
function(add_test_extrema TEST_NAME TEST_EXE PLTFILE TOL)

    setup_test()

    resolve_test_exe("${TEST_DIR}" "${TEST_EXE}" TEST_EXE)

    string(REPLACE ";" " " EXTREMA_ARGS "${ARGN}")
    set(test_command sh -c "${MPI_COMMANDS} ${TEST_EXE} ${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.i > ${TEST_NAME}.log && ${CMAKE_CURRENT_SOURCE_DIR}/check_extrema.sh ${FEXTREMA_EXE} ${CURRENT_TEST_BINARY_DIR}/${PLTFILE} ${TOL} ${EXTREMA_ARGS}")

    add_test(${TEST_NAME} ${test_command})
    set_tests_properties(${TEST_NAME}
        PROPERTIES
        TIMEOUT 5400
        PROCESSORS ${NP}
        WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
        LABELS "regression"
        ATTACHED_FILES_ON_FAIL "${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.log"
    )
endfunction(add_test_extrema)

# Stationary test -- compare with time 0
function(add_test_0 TEST_NAME TEST_EXE PLTFILE)
    setup_test()

    resolve_test_exe("${TEST_DIR}" "${TEST_EXE}" TEST_EXE)

    set(FCOMPARE_TOLERANCE "-r 1e-14 --abs_tol 1.0e-14")
    set(FCOMPARE_FLAGS "-a ${FCOMPARE_TOLERANCE}")
    set(test_command sh -c "${MPI_COMMANDS} ${TEST_EXE} ${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.i erf.input_sounding_file=${CURRENT_TEST_BINARY_DIR}/input_sounding ${RUNTIME_OPTIONS} > ${TEST_NAME}.log && ${FCOMPARE_EXE} ${FCOMPARE_FLAGS} ${CURRENT_TEST_BINARY_DIR}/plt00000 ${CURRENT_TEST_BINARY_DIR}/${PLTFILE}")

    add_test(${TEST_NAME} ${test_command})
    set_tests_properties(${TEST_NAME}
        PROPERTIES
        TIMEOUT 5400
        PROCESSORS ${NP}
        WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
        LABELS "regression"
        ATTACHED_FILES_ON_FAIL "${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.log"
    )
endfunction(add_test_0)

# Standard unit test
function(add_test_u TEST_NAME)
    setup_test()

    resolve_test_exe("${TEST_DIR}" "${TEST_EXE}" TEST_EXE)

    add_test(${TEST_NAME} sh -c "${MPI_COMMANDS} ${CMAKE_BINARY_DIR}/${amr_wind_unit_test_exe_name}")

    set_tests_properties(${TEST_NAME}
        PROPERTIES
        TIMEOUT 500
        PROCESSORS ${NP}
        WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
        LABELS "unit"
    )
endfunction(add_test_u)

#=============================================================================
# Unit tests
#=============================================================================
# add_test_u(unit_tests)

#=============================================================================
# Regression tests
#=============================================================================

add_test_r(DoublyPeriodic               "remora_exec" "plt00010")
add_test_r(DoublyPeriodic_bathy         "remora_exec" "plt00010")
add_test_r(Seamount                     "remora_exec" "plt00010")
add_test_r(Advection                    "remora_exec" "plt00010")
add_test_r(Advection_ML                 "remora_exec" "plt00010")
add_test_r(Upwelling                    "remora_exec" "plt00010")
add_test_r(Upwelling_GLS                "remora_exec" "plt00010")
add_test_r(Upwelling_NLEOS              "remora_exec" "plt00010")
add_test_r(Upwelling_qdrag              "remora_exec" "plt00010")
add_test_r(Upwelling_logdrag            "remora_exec" "plt00010")
add_test_r(Upwelling_Fennel             "remora_exec" "plt00010")
add_test_r(Channel_Test                 "remora_exec" "plt00010")
add_test_r(DoubleGyre                   "remora_exec" "plt00010")
add_test_r_hitol(BoundaryLayer          "remora_exec" "plt00010")
add_test_r(DogboneAnalytic              "remora_exec" "plt00010")
add_test_r(DogboneAnalytic_MLvel        "remora_exec" "plt_ml00010")
add_test_r(DogboneAnalytic_MLquad       "remora_exec" "plt_ml_quad00010")

#=============================================================================
# High-resolution initialization (remora.hires_grid_level / remora.hires_init_level)
#
# Bathymetry, grid metrics, or the initial state specified on a refined level and averaged
# down to level 0. The two transparency lanes reuse existing baselines, since a constant
# bathymetry must average down exactly; Seamount_hires is the lane that fails if the feature
# silently stops doing anything. hires_init_level is NetCDF-only, so it is covered by the
# developer lanes in Exec/GulfRefinementTest rather than here.
#=============================================================================

add_test_r_gold(Channel_Test_hires       "remora_exec" "plt00010"    Channel_Test)
add_test_r_gold(DogboneAnalytic_MLhires  "remora_exec" "plt_ml00010" DogboneAnalytic_MLvel)
add_test_r_differ(Seamount_hires         "remora_exec" "plt00010"    Seamount)
add_test_r_differ(Seamount_hires_r4      "remora_exec" "plt00010"    Seamount_hires)

# Six of the seven Fennel tracers are constants in the analytic profile, and a constant survives
# the average-down unchanged, so these are exact. tracer == 0 alongside them is what catches a
# Bio_comp = Tracer_comp + nscalar offset regression shifting biology into the dye slot.
add_test_extrema(Upwelling_Fennel_hires_init "remora_exec" "plt00000" 1e-12
                 tracer        0.0  0.0
                 NH4           0.1  0.1
                 chlorophyll   0.02 0.02
                 phytoplankton 0.08 0.08
                 zooplankton   0.06 0.06
                 LdetritusN    0.02 0.02
                 SdetritusN    0.04 0.04)

add_test_abort(Seamount_hires_init_abort      "remora_exec" "Cannot do high-resolution initialization for analytic initial conditions")
add_test_abort(Seamount_hires_grid_max_abort  "remora_exec" "hires_grid_level must be less than or equal to amr.max_level")
add_test_abort(Seamount_hires_init_max_abort  "remora_exec" "hires_init_level must be less than or equal to amr.max_level")
add_test_abort(Seamount_hires_grid_zero_abort "remora_exec" "hires_grid_level must be greater than 0")
add_test_abort(Seamount_hires_flat_abort      "remora_exec" "flat_bathymetry is incompatible with hires_grid_level")

#=============================================================================
# Performance tests
#=============================================================================

