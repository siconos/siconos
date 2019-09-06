# This file must be included in CMakeLists.txt when
# Siconos internal SuiteSparse is required (i.e. WITH_SYSTEM_SUITESPARSE=OFF).


# Compilation flags setup
include(CheckCCompilerFlag)
check_c_compiler_flag("-Wno-error=float-conversion" W_no_error_float_conversion_flag)
check_c_compiler_flag("-Wno-error=conversion" W_no_error_conversion_flag)
check_c_compiler_flag("-Wno-conversion" W_no_conversion_flag)
check_c_compiler_flag("-fpermissive" fpermissive_flag)

if(BUILD_AS_CPP)
  if(W_no_error_float_conversion_flag)
    list(APPEND _flags "-Wno-error=float-conversion")
  endif()
  if(W_no_error_conversion_flag)
    list(APPEND _flags "-Wno-error=conversion")
  endif()
  if(W_no_conversion_flag)
    list(APPEND _flags "-Wno-conversion")
  endif()
  if(fpermissive_flag)
    list(APPEND _flags "-fpermissive")
  endif()
else()
  if(W_no_error_conversion_flag)
    list(APPEND _flags "-Wno-error=conversion")
  endif()
  if(W_no_error_float_conversion_flag)
    list(APPEND _flags "-Wno-error=float-conversion")
  endif()
  if(W_no_conversion_flag)
    list(APPEND _flags "-Wno-conversion")
  endif()
endif()
list(REMOVE_DUPLICATES _flags)

# Set list of CXSparse files.
file(GLOB_RECURSE CXSPARSE_FILES ${CMAKE_CURRENT_SOURCE_DIR}  cxsparse_*.c)

# Add cxsparse sources to external build.
foreach(csfile IN LISTS CXSPARSE_FILES)
  target_sources(externals PRIVATE ${csfile})
  # set compilation flags
  foreach(flag IN LISTS _flags)
    set_source_files_properties(${csfile} PROPERTIES COMPILE_OPTIONS ${flag})
  endforeach()
endforeach()

