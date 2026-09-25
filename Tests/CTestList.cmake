
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
    # Validate the land/sea masks everywhere in the suite; it is off by default at run time.
    set(RUNTIME_OPTIONS "remora.check_mask_consistency=true")

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
    set(test_command sh -c "${MPI_COMMANDS} ${TEST_EXE} ${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.i ${RUNTIME_OPTIONS} > ${TEST_NAME}.log && ${FCOMPARE_EXE} ${FCOMPARE_FLAGS} ${FCOMPARE_GOLD_FILES_DIRECTORY}/${GOLD_NAME} ${CURRENT_TEST_BINARY_DIR}/${PLTFILE}")

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
    set(test_command sh -c "${MPI_COMMANDS} ${TEST_EXE} ${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.i ${RUNTIME_OPTIONS} > ${TEST_NAME}.log && ${FCOMPARE_EXE} ${FCOMPARE_FLAGS} ${PLOT_GOLD} ${CURRENT_TEST_BINARY_DIR}/${PLTFILE} && ! ${FCOMPARE_EXE} ${FCOMPARE_FLAGS} ${FCOMPARE_GOLD_FILES_DIRECTORY}/${OTHER_GOLD} ${CURRENT_TEST_BINARY_DIR}/${PLTFILE}")

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

# Two inputs describing the SAME configuration by different routes must agree. Neither run is a
# baseline, so a translation layer -- per-side against per-variable boundary specification -- is
# covered without a gold file blessed by the code under test, and neither route can drift alone.
# OTHER_INPUT.i sits beside TEST_NAME.i and must use a different plotfile prefix.
#
# Agreement alone cannot catch a value that both routes read wrongly in the SAME way: two runs
# that each ignore an input agree perfectly. Extra arguments, "<tol> <var> <min> <max> ...", add
# a check_extrema.sh assertion on TEST_NAME's plotfile, to pin the magnitude of something that
# must not go degenerate. Pick <tol> loose enough to read as an order-of-magnitude claim rather
# than a blessed digit string.
function(add_test_equiv TEST_NAME OTHER_INPUT TEST_EXE PLTFILE OTHER_PLTFILE)

    setup_test()

    resolve_test_exe("${TEST_DIR}" "${TEST_EXE}" TEST_EXE)

    set(FCOMPARE_TOLERANCE "-r 1e-11 --abs_tol 1.0e-11")
    set(FCOMPARE_FLAGS "-a ${FCOMPARE_TOLERANCE}")
    set(EXTREMA_CLAUSE "")
    if(NOT "${ARGN}" STREQUAL "")
        string(REPLACE ";" " " EXTREMA_ARGS "${ARGN}")
        set(EXTREMA_CLAUSE " && ${CMAKE_CURRENT_SOURCE_DIR}/check_extrema.sh ${FEXTREMA_EXE} ${CURRENT_TEST_BINARY_DIR}/${PLTFILE} ${EXTREMA_ARGS}")
    endif()
    set(test_command sh -c "${MPI_COMMANDS} ${TEST_EXE} ${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.i ${RUNTIME_OPTIONS} > ${TEST_NAME}.log 2>&1 && ${MPI_COMMANDS} ${TEST_EXE} ${CURRENT_TEST_BINARY_DIR}/${OTHER_INPUT}.i ${RUNTIME_OPTIONS} > ${OTHER_INPUT}.log 2>&1 && ${FCOMPARE_EXE} ${FCOMPARE_FLAGS} ${CURRENT_TEST_BINARY_DIR}/${OTHER_PLTFILE} ${CURRENT_TEST_BINARY_DIR}/${PLTFILE}${EXTREMA_CLAUSE}")

    add_test(${TEST_NAME} ${test_command})
    set_tests_properties(${TEST_NAME}
        PROPERTIES
        TIMEOUT 5400
        PROCESSORS ${NP}
        WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
        LABELS "regression"
        ATTACHED_FILES_ON_FAIL "${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.log;${CURRENT_TEST_BINARY_DIR}/${OTHER_INPUT}.log"
    )
endfunction(add_test_equiv)

# Test that a misconfigured input aborts, and aborts for the stated reason. The run must exit
# nonzero AND the log must carry the message, so a successful run, a different abort, and a
# segfault all fail. Deliberately not WILL_FAIL (which any nonzero exit satisfies) and not
# PASS_REGULAR_EXPRESSION (the model's stdout is redirected into the log, away from CTest).
# ABORT_SUBSTRING must be plain text. Note there is no ";" in the command: CMake would treat it
# as a list separator and hand sh -c two arguments, silently dropping everything after it.
function(add_test_abort TEST_NAME TEST_EXE ABORT_SUBSTRING)

    setup_test()

    resolve_test_exe("${TEST_DIR}" "${TEST_EXE}" TEST_EXE)

    set(test_command sh -c "! ${MPI_COMMANDS} ${TEST_EXE} ${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.i ${RUNTIME_OPTIONS} > ${TEST_NAME}.log 2>&1 && grep -q -- \"${ABORT_SUBSTRING}\" ${TEST_NAME}.log")

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

# Run the same input twice under different runtime options and require the two plotfiles to
# agree. For an invariant between two code paths -- neither run is a baseline, so this needs
# no gold data and cannot go stale against one.
function(add_test_r_selfcompare TEST_NAME TEST_EXE PLTFILE OPTIONS_A OPTIONS_B)

    setup_test()

    resolve_test_exe("${TEST_DIR}" "${TEST_EXE}" TEST_EXE)

    set(FCOMPARE_TOLERANCE "-r 1e-11 --abs_tol 1.0e-11")
    set(FCOMPARE_FLAGS "-a ${FCOMPARE_TOLERANCE}")
    set(test_command sh -c "mkdir -p runA runB && cd runA && ${MPI_COMMANDS} ${TEST_EXE} ${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.i ${OPTIONS_A} > ../${TEST_NAME}.log 2>&1 && cd ../runB && ${MPI_COMMANDS} ${TEST_EXE} ${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.i ${OPTIONS_B} >> ../${TEST_NAME}.log 2>&1 && cd .. && ${FCOMPARE_EXE} ${FCOMPARE_FLAGS} runA/${PLTFILE} runB/${PLTFILE}")

    add_test(${TEST_NAME} ${test_command})
    set_tests_properties(${TEST_NAME}
        PROPERTIES
        TIMEOUT 5400
        PROCESSORS ${NP}
        WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
        LABELS "regression"
        ATTACHED_FILES_ON_FAIL "${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.log"
    )
endfunction(add_test_r_selfcompare)

# Assert how far an integrated quantity drifts over a run. INPUT_NAME picks the input file, so
# several tests can share one case; OPTIONS go on the command line; the remaining arguments are
# <column> <bound> <below|above> triples handed to check_conservation.sh.
function(add_test_conservation TEST_NAME INPUT_NAME TEST_EXE OPTIONS)

    set(CURRENT_TEST_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${INPUT_NAME})
    set(CURRENT_TEST_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/test_files/${TEST_NAME})
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

    resolve_test_exe("${TEST_DIR}" "${TEST_EXE}" TEST_EXE)

    # sum_integrated_quantities returns immediately below verbosity 1, and its default six
    # digits cannot resolve the drifts asserted on here. The data log is opened for append, so
    # a stale one from an earlier run would supply the wrong first row.
    set(SUM_OPTS "remora.v=1 remora.sum_interval=1 remora.sum_precision=12 remora.data_log=cons.log")
    string(REPLACE ";" " " CHECKS "${ARGN}")

    set(test_command sh -c "rm -f cons.log && ${MPI_COMMANDS} ${TEST_EXE} ${CURRENT_TEST_BINARY_DIR}/${INPUT_NAME}.i ${SUM_OPTS} ${OPTIONS} > ${TEST_NAME}.log 2>&1 && ${CMAKE_CURRENT_SOURCE_DIR}/check_conservation.sh cons.log ${CHECKS}")

    add_test(${TEST_NAME} ${test_command})
    set_tests_properties(${TEST_NAME}
        PROPERTIES
        TIMEOUT 5400
        PROCESSORS ${NP}
        WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
        LABELS "regression"
        ATTACHED_FILES_ON_FAIL "${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.log"
    )
endfunction(add_test_conservation)

# Run must succeed AND its log must contain LOG_SUBSTRING. For a code path whose answers are
# not worth blessing into a gold file, but which must keep reaching the named behavior.
function(add_test_log TEST_NAME TEST_EXE LOG_SUBSTRING)

    setup_test()

    resolve_test_exe("${TEST_DIR}" "${TEST_EXE}" TEST_EXE)

    set(test_command sh -c "${MPI_COMMANDS} ${TEST_EXE} ${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.i > ${TEST_NAME}.log 2>&1 && grep -q -- \"${LOG_SUBSTRING}\" ${TEST_NAME}.log")

    add_test(${TEST_NAME} ${test_command})
    set_tests_properties(${TEST_NAME}
        PROPERTIES
        TIMEOUT 600
        PROCESSORS ${NP}
        WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
        LABELS "regression"
        ATTACHED_FILES_ON_FAIL "${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.log"
    )
endfunction(add_test_log)

# Assert variables' min and max in a plotfile against values known from outside REMORA -- a
# closed-form reference, or constants that an initial condition must reproduce. Unlike a gold
# file this says what the numbers should BE, so it also catches a baseline that was wrong when
# it was blessed. Pass any number of <var> <min> <max> triples after TOL; they are all checked
# against a single model run. See Tests/check_extrema.sh for the comparison itself.
function(add_test_extrema TEST_NAME TEST_EXE PLTFILE TOL)

    setup_test()

    resolve_test_exe("${TEST_DIR}" "${TEST_EXE}" TEST_EXE)

    string(REPLACE ";" " " EXTREMA_ARGS "${ARGN}")
    set(test_command sh -c "${MPI_COMMANDS} ${TEST_EXE} ${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.i ${RUNTIME_OPTIONS} > ${TEST_NAME}.log && ${CMAKE_CURRENT_SOURCE_DIR}/check_extrema.sh ${FEXTREMA_EXE} ${CURRENT_TEST_BINARY_DIR}/${PLTFILE} ${TOL} ${EXTREMA_ARGS}")

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

# Assert which refinement levels a run creates, and that it stays stationary while doing it.
# Grid structure is the only place a tagging bug shows: a criterion that has stopped tagging,
# or one that has started tagging the coastline, still writes a self-consistent plotfile, so
# fcompare alone cannot see it. MAX_LEVEL is the deepest level REMORA must build; 0 asserts
# that nothing was tagged. See Tests/check_max_level.sh for the assertion itself.
function(add_test_nlevels TEST_NAME TEST_EXE PLTFILE MAX_LEVEL)

    setup_test()

    resolve_test_exe("${TEST_DIR}" "${TEST_EXE}" TEST_EXE)

    set(FCOMPARE_TOLERANCE "-r 1e-14 --abs_tol 1.0e-14")
    set(FCOMPARE_FLAGS "-a ${FCOMPARE_TOLERANCE}")
    set(test_command sh -c "${MPI_COMMANDS} ${TEST_EXE} ${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.i ${RUNTIME_OPTIONS} > ${TEST_NAME}.log && ${CMAKE_CURRENT_SOURCE_DIR}/check_max_level.sh ${TEST_NAME}.log ${MAX_LEVEL} && ${FCOMPARE_EXE} ${FCOMPARE_FLAGS} ${CURRENT_TEST_BINARY_DIR}/plt00000 ${CURRENT_TEST_BINARY_DIR}/${PLTFILE}")

    add_test(${TEST_NAME} ${test_command})
    set_tests_properties(${TEST_NAME}
        PROPERTIES
        TIMEOUT 5400
        PROCESSORS ${NP}
        WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
        LABELS "regression"
        ATTACHED_FILES_ON_FAIL "${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.log"
    )
endfunction(add_test_nlevels)

# Assert how far a refinement level reaches, and that the run stays stationary. Where
# add_test_nlevels checks which levels exist, this checks where one of them is, which is the
# only thing that distinguishes "the region the user asked to refine was refined" from "part
# of it was". Pass the level, the direction (0/1/2 for x/y/z), and the physical extent it must
# span. See Tests/check_level_extent.sh.
function(add_test_extent TEST_NAME TEST_EXE PLTFILE LEVEL DIM LO HI)

    setup_test()

    resolve_test_exe("${TEST_DIR}" "${TEST_EXE}" TEST_EXE)

    set(FCOMPARE_TOLERANCE "-r 1e-14 --abs_tol 1.0e-14")
    set(FCOMPARE_FLAGS "-a ${FCOMPARE_TOLERANCE}")
    set(test_command sh -c "${MPI_COMMANDS} ${TEST_EXE} ${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.i ${RUNTIME_OPTIONS} > ${TEST_NAME}.log && ${CMAKE_CURRENT_SOURCE_DIR}/check_level_extent.sh ${CURRENT_TEST_BINARY_DIR}/${PLTFILE} ${LEVEL} ${DIM} ${LO} ${HI} && ${FCOMPARE_EXE} ${FCOMPARE_FLAGS} ${CURRENT_TEST_BINARY_DIR}/plt00000 ${CURRENT_TEST_BINARY_DIR}/${PLTFILE}")

    add_test(${TEST_NAME} ${test_command})
    set_tests_properties(${TEST_NAME}
        PROPERTIES
        TIMEOUT 5400
        PROCESSORS ${NP}
        WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
        LABELS "regression"
        ATTACHED_FILES_ON_FAIL "${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.log"
    )
endfunction(add_test_extent)

# Assert how many cells a refinement level covers. The third grid assertion, for a change that
# removes cells from inside a refined region while leaving its bounding box alone -- which
# neither add_test_nlevels nor add_test_extent can see. Summed over boxes, so it does not
# depend on how the region is chopped for load balance. No stationary comparison: use this for
# a case that is actually evolving. See Tests/check_level_cells.sh.
function(add_test_cells TEST_NAME TEST_EXE PLTFILE LEVEL NCELLS)

    setup_test()

    resolve_test_exe("${TEST_DIR}" "${TEST_EXE}" TEST_EXE)

    set(test_command sh -c "${MPI_COMMANDS} ${TEST_EXE} ${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.i ${RUNTIME_OPTIONS} > ${TEST_NAME}.log && ${CMAKE_CURRENT_SOURCE_DIR}/check_level_cells.sh ${CURRENT_TEST_BINARY_DIR}/${PLTFILE} ${LEVEL} ${NCELLS}")

    add_test(${TEST_NAME} ${test_command})
    set_tests_properties(${TEST_NAME}
        PROPERTIES
        TIMEOUT 5400
        PROCESSORS ${NP}
        WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
        LABELS "regression"
        ATTACHED_FILES_ON_FAIL "${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.log"
    )
endfunction(add_test_cells)

# Stationary test -- compare with time 0
function(add_test_0 TEST_NAME TEST_EXE PLTFILE)
    setup_test()

    resolve_test_exe("${TEST_DIR}" "${TEST_EXE}" TEST_EXE)

    set(FCOMPARE_TOLERANCE "-r 1e-14 --abs_tol 1.0e-14")
    set(FCOMPARE_FLAGS "-a ${FCOMPARE_TOLERANCE}")
    set(test_command sh -c "${MPI_COMMANDS} ${TEST_EXE} ${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.i ${RUNTIME_OPTIONS} > ${TEST_NAME}.log && ${FCOMPARE_EXE} ${FCOMPARE_FLAGS} ${CURRENT_TEST_BINARY_DIR}/plt00000 ${CURRENT_TEST_BINARY_DIR}/${PLTFILE}")

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
# A refined patch whose low x edge lies on the periodic seam. The only lane whose fine level
# touches a periodic boundary, so the only one that holds the periodic-boundary handling in
# BuildMask, fill_ghost_kcomps, FillCoarsePatchMap and the ubar/vbar average-down in place.
# Also the only multi-level lane with Coriolis on. Its gold is from one rank, one box per
# level; under MPI the level-1 box splits at j = 30, so it also holds the interface ghost
# refresh after set_2d_cf_bcs, without which the answer depended on the box layout.
add_test_r(Channel_Test_ML_seam         "remora_exec" "plt00010")

#=============================================================================
# Time subcycling on refined levels (remora.do_substep)
#
# The default path's answers are pinned by Advection_ML and the two DogboneAnalytic_ML gold
# lanes. These assert behaviour those cannot: that subcycling actually engages, and that the
# two drivers still agree where they should. Here, level 1 receives dt[0]/2 (fixed_dt = 100,
# ref_ratio = 2).
#=============================================================================
add_test_log(Advection_ML_subcycle      "remora_exec" "with dt = 50")

# The load-bearing one: at a timestep ratio of 1 the recursive driver must reproduce
# timeStepML, separating a broken driver from the answer changes subcycling legitimately
# makes. With DogboneAnalytic_ML_conservation_lockstep it is one of the two things still
# holding lockstep in place.
# do_reflux is off in the subcycled run because it is a correction the lockstep driver does
# not apply at all, so leaving it on would compare a feature rather than the drivers. It
# moves the tracer by 3e-4 here, well clear of the tolerance.
add_test_r_selfcompare(Advection_ML_subcycle_identity "remora_exec" "plt00020"
                       "remora.do_substep=0"
                       "remora.do_substep=1 remora.dt_ref_ratio=1 remora.do_reflux=0")

# Advection has a flat bottom, so D matches across the interface and set_2d_cf_bcs reduces to
# the interpolation it replaces. This lane has varying bathymetry and a refinement ratio of 3,
# so the mass-flux form is actually exercised.
add_test_log(DogboneAnalytic_ML_subcycle "remora_exec" "3 x 3          3          60       0.6667")

#=============================================================================
# Conservation
#
# Advection is doubly periodic and DogboneAnalytic is closed by slipwalls, so in both nothing
# can leave the domain and the totals have to hold. Bounds come from measurement, not from
# taste; the numbers each one is separating are in the comments.
#=============================================================================

# Tracer mass. Refluxing takes the drift from 7.8e-5 to below what 12 digits can resolve.
add_test_conservation(Advection_ML_conservation Advection_ML_subcycle "remora_exec"
                      "remora.max_step=20 remora.do_reflux=1 remora.reflux_clamp=0"
                      tracer 1e-10 below)

# The control, and the reason the lane above means anything: without refluxing the same run
# must drift. If this ever passes by conserving, the case has stopped exercising the
# correction -- no interface, no gradient across it, or a no-op -- and its partner above is
# proving nothing.
add_test_conservation(Advection_ML_conservation_control Advection_ML_subcycle "remora_exec"
                      "remora.max_step=20 remora.do_reflux=0"
                      tracer 1e-6 above)

# Three levels, where the correction at the 1/2 interface has to be accumulated over several
# steps of level 1 before it is applied. A two-level case refluxes once per step of the only
# coarse level there is, so it passes whether or not that accumulation is right.
add_test_conservation(Advection_3L_conservation Advection_3L_conservation "remora_exec"
                      "remora.max_step=20 remora.do_reflux=1 remora.reflux_clamp=0"
                      tracer 1e-10 below)

# Its control, for the same reason as above.
add_test_conservation(Advection_3L_conservation_control Advection_3L_conservation "remora_exec"
                      "remora.max_step=20 remora.do_reflux=0"
                      tracer 1e-6 above)

# The floor. Both totals are exact on one level, so the AMR bounds are measured against
# roundoff rather than against an unknown scheme error.
add_test_conservation(Advection_conservation_baseline Advection_ML_subcycle "remora_exec"
                      "remora.max_step=20 amr.max_level=0"
                      tracer 1e-12 below volume 1e-12 below)

# Volume is the integral of the free surface, so what costs it here is how that surface is
# written in the fine ghost band. REMORA writes it as ROMS's put_refine2d does -- every
# leapfrog record plus Zt_avg1 (cf_fill_all_kcomp), interpolated onto the child's sub-time
# (cf_time_interp_zeta) -- and neither preserves volume: one at a time they take the drift
# from 2.0e-09 to 2.9e-08 and 1.4e-08 (measured; no lane turns off one alone). That is the
# trade ROMS makes. The default path measures 2.9e-08, against 2.0e-09 for the lane below and
# 2.0e-6 for the lockstep driver, which interpolates ubar instead of imposing the flux. Single
# level is exact on this case too, measured; the lane asserting that is the Advection one.
add_test_conservation(DogboneAnalytic_ML_conservation DogboneAnalytic_ML_subcycle "remora_exec"
                      "remora.max_step=20"
                      volume 1e-7 below)

# The conservative path, pinned where the bound above used to sit. Switching off only the two
# ghost-band knobs keeps ROMS's restriction stencil and flux distribution and brings the drift
# back to 2.0e-09, where the transport side of the interface sits on its own: of that,
# cf_avgdown_stencil's nine-point mean costs about 3e-10 and cf_flux_pc is exactly neutral
# (both measured; no lane isolates either). If this lane drifts, the conservative path has
# regressed whatever the default path is doing.
add_test_conservation(DogboneAnalytic_ML_conservation_cons DogboneAnalytic_ML_subcycle "remora_exec"
                      "remora.max_step=20 remora.cf_fill_all_kcomp=0 remora.cf_time_interp_zeta=0"
                      volume 1e-8 below)

# The must-drift control for both volume lanes above, and the assertion behind the lockstep
# figure they are quoted against. Lockstep interpolates the parent's ubar at the interface
# instead of imposing its fast-time-averaged flux -- its shared barotropic loop leaves DU_avg2
# still accumulating when a child needs it -- and volume drifts 2.0e-6, two orders past the
# 1e-7 the default path is held to. If this lane ever conserves, the case has stopped
# separating the two drivers and neither bound above is measuring the coarse-fine treatment.
add_test_conservation(DogboneAnalytic_ML_conservation_lockstep DogboneAnalytic_ML_subcycle "remora_exec"
                      "remora.max_step=20 remora.do_substep=0"
                      volume 1e-7 above)

# The assumption underneath all of the above: the fine cell edges have to sum to the coarse
# edge, or the mass flux imposed at the interface cannot be conservative whatever else is
# right. check_cf_metrics aborts past remora.check_cf_tol, so reaching the printed line is the
# assertion. It measures 0 on this ratio-2 analytic grid and 1.4e-16 on Dogbone's ratio 3, but
# is not guaranteed on the NetCDF path, where a finer level interpolates its metrics from the
# parent's and scales them by the refinement ratio. All three lanes set the tolerance to 1e-14,
# within a factor of 40 of the largest measurement, rather than the 1e-12 default.
add_test_log(Advection_ML_cf_metrics "remora_exec" "CF edge tiling")

# The ratio-3 grid the DogboneAnalytic_ML golds and both volume lanes above are measured on,
# so the identity those numbers rest on is asserted where they are taken rather than only at
# ratio 2: three fine faces to a coarse one, 1.4e-16.
add_test_log(DogboneAnalytic_ML_cf_metrics "remora_exec" "CF edge tiling")

# The same check where it is not trivially satisfied: every other multi-level case has uniform
# pm and pn, so their fine edges sum to the coarse edge for a reason particular to them.
add_test_log(BoundaryLayer_ML_cf_metrics "remora_exec" "CF edge tiling")

# amr.do_substep is the original spelling and has to keep working. The warning is the
# observable proof the fallback was read rather than silently ignored, and setting both
# spellings is an error rather than a silent precedence rule.
add_test_log(Advection_ML_do_substep_alias "remora_exec" "amr.do_substep is deprecated")
add_test_abort(Advection_ML_do_substep_both_abort "remora_exec" "and amr.do_substep are both")

#=============================================================================
# High-resolution initialization (remora.hires_grid_level / remora.hires_init_level)
#
# Bathymetry, grid metrics, or the initial state specified on a refined level and averaged
# down to level 0. The two transparency lanes reuse existing baselines, since a constant
# bathymetry must average down exactly; Seamount_hires is the lane that fails if the feature
# silently stops doing anything. hires_init_level is covered here for analytic initial
# conditions, and its NetCDF reader by the developer lanes in Exec/GulfRefinementTest.
#=============================================================================

add_test_r_gold(Channel_Test_hires       "remora_exec" "plt00010"    Channel_Test)
add_test_r_gold(DogboneAnalytic_MLhires  "remora_exec" "plt_ml00010" DogboneAnalytic_MLvel)

# The only tests with a partially-masked coarse cell. Their exact solution is rest, so
# they need no gold file: plt00010 must equal plt00000. A plain arithmetic average-down
# lets the zeroed fine land cells drag those coarse cells off their initial value, which
# breaks stationarity by ~1e-2 in salt and velocity.
# coastbox spans the full width of the domain, 0 to 750 in y, and the dogbone's coasts run
# through it, so this asserts the thing a static box is for: the whole box is refined, land
# included. REMORA used to clip it back to the water, giving y in [250,500].
add_test_extent(DogboneAnalytic_MLmask  "remora_exec" "plt00010" 1 1 0.0 750.0)
# Same assertion on the ref-ratio-2 twin, so that path of the same fix is covered too: its
# grids move from 3 clipped boxes to 2 full-box grids, which add_test_0 alone cannot see.
add_test_extent(DogboneAnalytic_MLmask_rr2 "remora_exec" "plt00010" 1 1 0.0 750.0)
# A pair: the same gradient criterion keyed on a physical field, which must not see the coast,
# and on the mask, which must. Either one alone passes for the wrong reason.
add_test_nlevels(DogboneAnalytic_MLcoastskip "remora_exec" "plt00010" 0)
add_test_nlevels(DogboneAnalytic_MLcoasttag  "remora_exec" "plt00010" 1)
# MLcoastskip covers a cell-centered field, where the cell's own mask is the whole story.
# This covers a face-staggered one, where it is not: 12960 cells are refined if x_velocity is
# guarded by mskr(i,j) alone instead of by both cells sharing its face.
add_test_cells(DogboneAnalytic_MLdryface     "remora_exec" "plt00010" 1 6480)
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

# Analytic hires_init_level. Seamount's temperature is convex in z and z_r varies within a
# coarse cell, so averaging down differs from evaluating at level 0 (max temp 17.475031408 with
# the flag off). Expected extrema come from Tests/reference/hires_init_analytic_reference.py,
# which matches the whole level-0 field to 1e-14. One lane per bathymetry source at the init
# level: analytic there, analytic there above coarser grid data, averaged down from finer.
add_test_extrema(Seamount_hires_init             "remora_exec" "plt00000" 1e-8
                 temp 10.077862333 17.474890098
                 salt 32.0 32.0)
add_test_extrema(Seamount_hires_init_above_grid  "remora_exec" "plt00000" 1e-8
                 temp 10.077862333 17.474877547
                 salt 32.0 32.0)
add_test_extrema(Seamount_hires_init_below_grid  "remora_exec" "plt00000" 1e-8
                 temp 10.077862333 17.474875997
                 salt 32.0 32.0)

# Rest must stay rest when the initial state is also averaged down through partial masks.
add_test_0(DogboneAnalytic_MLmask_hires_init    "remora_exec" "plt00010")

# Biology evaluated on the hires level and averaged down: the constants survive exactly, and an
# unfilled component would come out 0. The max temp is a snapshot; it separates this lane from
# Upwelling_Fennel_hires_init (21.935181561), which does not set hires_init_level.
add_test_extrema(Upwelling_Fennel_hires_init_ic "remora_exec" "plt00000" 1e-8
                 temp          14.520483402 21.935189849
                 tracer        0.0  0.0
                 NH4           0.1  0.1
                 chlorophyll   0.02 0.02
                 phytoplankton 0.08 0.08
                 zooplankton   0.06 0.06
                 LdetritusN    0.02 0.02
                 SdetritusN    0.04 0.04)

add_test_abort(Seamount_hires_init_gridscale_abort "remora_exec" "requires remora.grid_scale_type = constant")
add_test_abort(Seamount_hires_grid_max_abort  "remora_exec" "hires_grid_level must be less than or equal to amr.max_level")
add_test_abort(Seamount_hires_init_max_abort  "remora_exec" "hires_init_level must be less than or equal to amr.max_level")
add_test_abort(Seamount_hires_grid_zero_abort "remora_exec" "hires_grid_level must be greater than 0")

#=============================================================================
# Boundary conditions
#
# BC_per_variable states four conditions as one West South East North list and must land them on
# the right faces; BC_per_side names each face outright and cannot get the order wrong.
#
# The y pairing is pinned by the VELOCITY conditions: for cell-centered tracers, and for zeta and
# tke, noslipwall and slipwall both map to foextrap, so South and North are indistinguishable in
# those fields. Only xvel/yvel and ubar/vbar tell them apart -- no_slip_wall gives ext_dir on both
# components, slip_wall foextrap tangential and ext_dir normal. Swapping in conditions that look
# equally distinct but agree for the velocities would silently blind this lane.
#
# The dye is zero in the initial condition and enters only through the western inflow value, so it
# is nonzero only if that value was read. Agreement alone cannot assert that -- two runs that both
# ignored the value would agree at zero -- so the extrema clause pins its magnitude. The tolerance
# is deliberately loose: the claim is "order 0.1, and certainly not zero", not a blessed digit
# string. Entry is by horizontal diffusion; see BC_per_variable.i on why it is not advective.
#=============================================================================

add_test_equiv(BC_per_variable BC_per_side "remora_exec" "plt00010" "plt_side00010"
               0.05 tracer 0.0 0.127)

add_test_abort(BC_inflow_no_value "remora_exec" "needs an inflow value")

#=============================================================================
# Barotropic substep count
#
# remora.ndtfast is the number of barotropic steps per baroclinic step. Advance and
# timeStepML form the fast step as dt / ndtfast and set_weights sizes the barotropic filter
# with it, none of which is guarded at the point of use, so it has to be positive by the time
# parsing finishes. It has no usable default: it was once inferred from remora.fixed_dt /
# remora.fixed_fast_dt, which left it at zero on every path that did not set both -- a
# CFL-driven run cannot set them, since dt is not known until run time. Zero divides by zero
# and sizes the weight vectors to one element, which advance_2d then reads past the end of;
# amrex::Vector does not bounds-check that. So the unset case must abort, and does.
#
# The other three lanes cover the input names around it. fixed_fast_dt is gone and must be
# rejected rather than ignored: amrex does not abort on unused inputs, so a stale input file
# naming it would otherwise run with whatever substep count it did not ask for. The
# fixed_ndtfast_ratio alias must still deliver the count it names -- agreement here cannot
# pass for the wrong reason, because a run that ignored the alias would be left at zero and
# abort rather than agree -- and naming the count under both spellings at once is an error,
# since silently preferring one leaves the input file reading as though it asked for the other.
#=============================================================================

add_test_equiv(Channel_Test_ndtfast_alias Channel_Test_ndtfast_named "remora_exec" "plt00010" "plt_named00010")

add_test_abort(Channel_Test_ndtfast_unset_abort  "remora_exec" "remora.ndtfast must be a positive integer")
add_test_abort(Channel_Test_fixed_fast_dt_abort  "remora_exec" "remora.fixed_fast_dt has been removed")
add_test_abort(Channel_Test_ndtfast_both_abort   "remora_exec" "and remora.fixed_ndtfast_ratio are both")

#=============================================================================
# Performance tests
#=============================================================================

