set(remora_exe_name advection)

add_executable(${remora_exe_name} "")
target_sources(${remora_exe_name}
   PRIVATE
     prob.H
     prob.cpp
)

target_include_directories(${remora_exe_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

include(${CMAKE_SOURCE_DIR}/CMake/BuildREMORAExe.cmake)
build_remora_exe(${remora_exe_name})
