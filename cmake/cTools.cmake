# =====================================
# Some macros to check C compiler,
# check/add some flags, get version ...
# =====================================

include(CheckCCompilerFlag)

# --- Append some options (list OPT) to cxx compiler flags. ---
# Update :
# - CMAKE_C_FLAGS
# -
# parameters :
# OPT : the flags
# extra parameters : a list of compilers ids
# (see http://www.cmake.org/cmake/help/v3.3/variable/CMAKE_LANG_COMPILER_ID.html#variable:CMAKE_<LANG>_COMPILER_ID)
# If extra parameter is set, the option will be applied
# only for those compilers, if the option is accepted.
# 
macro(ADD_C_OPTIONS OPT)
  
  STRING(REGEX REPLACE " " "" OPT_SANE "${OPT}")
  CHECK_C_COMPILER_FLAG("${OPT} ${_EXTRA_WARNING_FLAGS}" C_HAVE_${OPT_SANE})
  set(_compilers ${ARGN})
  IF(_compilers)
    SET(ADD_OPTION FALSE)
    FOREACH(_compiler ${_compilers})
      IF (${CMAKE_C_COMPILER_ID} MATCHES ${_compiler})
	MESSAGE(STATUS "Adding option for compiler ${_compiler}")
	SET(ADD_OPTION TRUE)
      ENDIF()
    ENDFOREACH()
  ELSE(_compilers)
    SET(ADD_OPTION TRUE)
  ENDIF(_compilers)
  
  IF(ADD_OPTION AND C_HAVE_${OPT_SANE})
    APPEND_C_FLAGS("${OPT}")
  ENDIF(ADD_OPTION AND C_HAVE_${OPT_SANE})
endmacro()

# Based on the Qt 5 processor detection code, so should be very accurate
# https://qt.gitorious.org/qt/qtbase/blobs/master/src/corelib/global/qprocessordetection.h
# Currently handles arm (v5, v6, v7), x86 (32/64), ia64, and ppc (32/64)

# Regarding POWER/PowerPC, just as is noted in the Qt source,
# "There are many more known variants/revisions that we do not handle/detect."

set(c_detect_code "
#if __STDC_VERSION__ >= 201112L
#error STDC_VERSION 2011012L
#else
#error STDC_VERSION 000000L
#endif                                                                                                                                                                                                                                      ")

# Set ppc_support to TRUE before including this file or ppc and ppc64
# will be treated as invalid architectures since they are no longer supported by Apple
function(detect_c_version output_var)
  file(WRITE "${CMAKE_BINARY_DIR}/cversion.c" "${c_detect_code}")
  
  enable_language(C)
  try_run(
    run_result_unused
    compile_result_unused
    "${CMAKE_BINARY_DIR}"
    "${CMAKE_BINARY_DIR}/cversion.c"
    COMPILE_OUTPUT_VARIABLE CVERSION
    )
  
  # Parse the architecture name from the compiler output
  string(REGEX MATCH "STDC_VERSION ([a-zA-Z0-9_]+)" CVERSION "${CVERSION}")
  
  # Get rid of the value marker leaving just the architecture name
  string(REPLACE "STDC_VERSION " "" CVERSION "${CVERSION}")
  
  set(${output_var} "${CVERSION}" PARENT_SCOPE)
endfunction()
