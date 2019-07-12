# Switch between system or local suitesparse,
# depending on the chosen user option.

# Add sources dir for suite parse (local)
list(APPEND ${COMPONENT}_DIRS SuiteSparse/CXSparse)
set_source_files_properties(SuiteSparse/CSparse/csparse.c PROPERTIES COMPILE_OPTIONS -Wno-unused)

include(CheckCCompilerFlag)
check_c_compiler_flag("-Wno-error=float-conversion" W_no_error_float_conversion_flag)
check_c_compiler_flag("-Wno-error=conversion" W_no_error_conversion_flag)
check_c_compiler_flag("-Wno-conversion" W_no_conversion_flag)
check_c_compiler_flag("-fpermissive" fpermissive_flag)

file(GLOB_RECURSE CXSPARSE_FILES ${CMAKE_CURRENT_SOURCE_DIR}  cxsparse_*.c)
if(BUILD_AS_CPP)
  set(_flags)
  if(W_no_error_float_conversion_flag)
    set(_flags "-Wno-error=float-conversion")
  endif()
  if(W_no_error_conversion_flag)
    set(_flags "${_flags} -Wno-error=conversion")
  endif()
  if(W_no_conversion_flag)
    set(_flags "${_flags} -Wno-conversion")
  endif()
  if(fpermissive_flag)
    set(_flags "${_flags} -fpermissive")
  endif()
  if(_flags)
    set_source_files_properties(${CXSPARSE_FILES} PROPERTIES COMPILE_OPTIONS ${_flags})
  endif()
else()
  set(_flags)
  if(W_no_error_conversion_flag)
    set(_flags "-Wno-error=conversion")
  endif()
  if(W_no_error_float_conversion_flag)
    set(_flags "${_flags} -Wno-error=float-conversion")
  endif()
  if(W_no_conversion_flag)
    set(_flags "${_flags} -Wno-conversion")
  endif()
  if(_flags)
    set_source_files_properties(${CXSPARSE_FILES} PROPERTIES COMPILE_OPTIONS ${_flags})
  endif()
endif()