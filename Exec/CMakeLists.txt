set(remora_lib_name remora_srclib)

add_library(${remora_lib_name} OBJECT)
if (BUILD_SHARED_LIBS)
  set_target_properties(${remora_lib_name} PROPERTIES POSITION_INDEPENDENT_CODE ON)
endif()

include(${PROJECT_SOURCE_DIR}/CMake/BuildREMORAExe.cmake)

build_remora_lib(${remora_lib_name})

if(NOT REMORA_BUILD_LIBRARY_ONLY)
  set(remora_exe_name remora_exec)
  add_executable(${remora_exe_name} "")
  target_sources(${remora_exe_name}
     PRIVATE
       ${PROJECT_SOURCE_DIR}/Exec/REMORA_Prob.cpp
  )
  target_include_directories(${remora_exe_name} PRIVATE ${PROJECT_SOURCE_DIR}/Exec)
  build_remora_exe(${remora_exe_name})

endif() # NOT REMORA_BUILD_LIBRARY_ONLY
