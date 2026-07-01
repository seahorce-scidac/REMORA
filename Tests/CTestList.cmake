
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

# Startup contract test: run a short case and assert expected log strings exist.
function(add_test_log_contains TEST_NAME TEST_EXE RUN_ARGS)
    setup_test()

    set(TEST_EXE ${CMAKE_BINARY_DIR}/Exec/${TEST_EXE})
    set(test_command sh -c "${MPI_COMMANDS} ${TEST_EXE} ${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.i ${RUN_ARGS} > ${TEST_NAME}.log")

    add_test(${TEST_NAME} ${test_command})
    set_tests_properties(${TEST_NAME}
        PROPERTIES
        TIMEOUT 1200
        PROCESSORS ${NP}
        WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
        LABELS "startup_contract"
        ATTACHED_FILES_ON_FAIL "${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.log"
    )

    foreach(expected_line IN LISTS ARGN)
        string(REGEX REPLACE "[^A-Za-z0-9_]" "_" expected_tag "${expected_line}")
        set(check_name "${TEST_NAME}_check_${expected_tag}")
        add_test(${check_name} sh -c "grep -Fq \"${expected_line}\" ${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.log")
        set_tests_properties(${check_name}
            PROPERTIES
            TIMEOUT 120
            PROCESSORS 1
            WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
            LABELS "startup_contract"
            ATTACHED_FILES_ON_FAIL "${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.log"
        )
        set_tests_properties(${check_name} PROPERTIES DEPENDS ${TEST_NAME})
    endforeach()
endfunction(add_test_log_contains)

# Startup contract failure test: command must fail and log must contain expected strings.
function(add_test_log_contains_fail TEST_NAME TEST_EXE RUN_ARGS)
    setup_test()

    set(TEST_EXE ${CMAKE_BINARY_DIR}/Exec/${TEST_EXE})
    set(test_command sh -c "${MPI_COMMANDS} ${TEST_EXE} ${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.i ${RUN_ARGS} > ${TEST_NAME}.log 2>&1")

    add_test(${TEST_NAME} ${test_command})
    set_tests_properties(${TEST_NAME}
        PROPERTIES
        WILL_FAIL TRUE
        TIMEOUT 1200
        PROCESSORS ${NP}
        WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
        LABELS "startup_contract"
        ATTACHED_FILES_ON_FAIL "${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.log"
    )

    foreach(expected_line IN LISTS ARGN)
        string(REGEX REPLACE "[^A-Za-z0-9_]" "_" expected_tag "${expected_line}")
        set(check_name "${TEST_NAME}_check_${expected_tag}")
        add_test(${check_name} sh -c "grep -Fq \"${expected_line}\" ${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.log")
        set_tests_properties(${check_name}
            PROPERTIES
            TIMEOUT 120
            PROCESSORS 1
            WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
            LABELS "startup_contract"
            ATTACHED_FILES_ON_FAIL "${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.log"
        )
        set_tests_properties(${check_name} PROPERTIES DEPENDS ${TEST_NAME})
    endforeach()
endfunction(add_test_log_contains_fail)

#=============================================================================
# Unit tests
#=============================================================================
# add_test_u(unit_tests)

#=============================================================================
# Regression tests
#=============================================================================

    add_test_log_contains(BiologyStartup_ECB "DoublyPeriodic/doublyperiodic"
            "remora.max_step=0 remora.sum_interval=1 remora.v=1 remora.biology_enabled=true remora.biology_model=ecb remora.biology_require_tracers_strict=true remora.nscalar=1"
            "[REMORA] Expanding nscalar from 1 to 9 to satisfy ecb packed tracer requirements through DON"
            "[REMORA] Biology tracer registry entries: 12"
            "[REMORA] Biology startup status: enabled=true, model=ecb, strict_tracer_check=true"
            "biology_scaffold_required_total=9")

    add_test_log_contains(BiologyStartup_Fennel "DoublyPeriodic/doublyperiodic"
            "remora.max_step=0 remora.sum_interval=1 remora.v=1 remora.biology_enabled=true remora.biology_model=fennel remora.biology_require_tracers_strict=true remora.nscalar=1"
            "[REMORA] Expanding nscalar from 1 to 9 to satisfy fennel packed tracer requirements through DON"
            "[REMORA] Biology tracer registry entries: 12"
            "[REMORA] Biology startup status: enabled=true, model=fennel, strict_tracer_check=true"
            "biology_scaffold_required_total=9")

      add_test_log_contains(BiologyStartup_Ecosim "DoublyPeriodic/doublyperiodic"
          "remora.max_step=0 remora.sum_interval=1 remora.v=1 remora.biology_enabled=true remora.biology_model=ecosim remora.biology_require_tracers_strict=true remora.nscalar=1"
          "[REMORA] Expanding nscalar from 1 to 9 to satisfy ecosim packed tracer requirements through DON"
          "[REMORA] Biology tracer registry entries: 12"
          "[REMORA] Biology startup status: enabled=true, model=ecosim, strict_tracer_check=true"
          "biology_scaffold_required_total=9")

      add_test_log_contains(BiologyStartup_NPZDFranks "DoublyPeriodic/doublyperiodic"
          "remora.max_step=0 remora.sum_interval=1 remora.v=1 remora.biology_enabled=true remora.biology_model=npzd_franks remora.biology_require_tracers_strict=true remora.nscalar=1"
          "[REMORA] Biology tracer registry entries: 0"
          "[REMORA] Biology startup status: enabled=true, model=npzd_franks, strict_tracer_check=true"
          "biology_scaffold_required_total=4"
          "biology_scaffold_missing_required_flag=1")

      add_test_log_contains(BiologyStartup_NPZDFranksAlias "DoublyPeriodic/doublyperiodic"
          "remora.max_step=0 remora.sum_interval=1 remora.v=1 remora.biology_enabled=true remora.biology_model=npzd-franks remora.biology_require_tracers_strict=true remora.nscalar=1"
          "[REMORA] Biology tracer registry entries: 0"
          "[REMORA] Biology startup status: enabled=true, model=npzd-franks, strict_tracer_check=true"
          "[REMORA] Creating NPZD-Franks biology model plugin"
          "biology_scaffold_required_total=4"
          "biology_scaffold_missing_required_flag=1")

      add_test_log_contains_fail(BiologyStartup_InvalidNscalar "DoublyPeriodic/doublyperiodic"
          "remora.max_step=0 remora.v=1 remora.nscalar=0"
          "remora.nscalar must be at least 1")
endif()

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
add_test_r(Channel_Test                 "remora_exec" "plt00010")
add_test_r(DoubleGyre                   "remora_exec" "plt00010")
add_test_r_hitol(BoundaryLayer          "remora_exec" "plt00010")
add_test_r(DogboneAnalytic              "remora_exec" "plt00010")
add_test_r(DogboneAnalytic_MLvel        "remora_exec" "plt_ml00010")
add_test_r(DogboneAnalytic_MLquad       "remora_exec" "plt_ml_quad00010")

#=============================================================================
# Performance tests
#=============================================================================

