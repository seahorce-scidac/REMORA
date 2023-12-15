set(romsx_exe_name doublyperiodic)

add_executable(${romsx_exe_name} "")
target_sources(${romsx_exe_name}
   PRIVATE
     prob.H
     prob.cpp
)

target_include_directories(${romsx_exe_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

include(${CMAKE_SOURCE_DIR}/CMake/BuildROMSXExe.cmake)
build_romsx_exe(${romsx_exe_name})
