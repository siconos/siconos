# =====================================
# Some macros to check CXX compiler,
# check/add some flags, get version ...
# =====================================

include(TestCXXAcceptsFlag)

# --- Append some options (list OPT) to cxx compiler flags. ---
# Update :
# - CMAKE_CXX_FLAGS
# -
# parameters :
# OPT : the flags
# extra parameters : a list of compilers ids
# (see http://www.cmake.org/cmake/help/v3.3/variable/CMAKE_LANG_COMPILER_ID.html#variable:CMAKE_<LANG>_COMPILER_ID)
# If extra parameter is set, the option will be applied
# only for those compilers, if the option is accepted.
# 
macro(ADD_CXX_OPTIONS OPT)
  STRING(REGEX REPLACE " " "" OPT_SANE "${OPT}")
  # Check if option exists for the current compiler ...
  CHECK_CXX_ACCEPTS_FLAG("${OPT}" CXX_HAVE_${OPT_SANE})

  # deal with extra arguments (others than OPT list)
  set(_compilers ${ARGN})
  IF(_compilers)
    SET(ADD_OPTION FALSE)
    FOREACH(_compiler ${_compilers})
      IF (${CMAKE_CXX_COMPILER_ID} MATCHES ${_compiler})
	SET(ADD_OPTION TRUE)
      ENDIF()
    ENDFOREACH()
  ELSE(_compilers)
    SET(ADD_OPTION TRUE)
  ENDIF(_compilers)
  IF(ADD_OPTION AND CXX_HAVE_${OPT_SANE})
    APPEND_CXX_FLAGS("${OPT}")
  ENDIF(ADD_OPTION AND CXX_HAVE_${OPT_SANE})
endmacro(ADD_CXX_OPTIONS)

# Based on the Qt 5 processor detection code, so should be very accurate
# https://qt.gitorious.org/qt/qtbase/blobs/master/src/corelib/global/qprocessordetection.h
# Currently handles arm (v5, v6, v7), x86 (32/64), ia64, and ppc (32/64)

# Regarding POWER/PowerPC, just as is noted in the Qt source,
# "There are many more known variants/revisions that we do not handle/detect."

set(cxx_detect_code "
#if __cplusplus >= 201103L
#error CXXVERSION 201103L
#elif __cplusplus >= 199711L
#error CXXVERSION 199711L
#else
#error CXXVERSION 000000L
#endif
")

# --------------- Return Cxx compiler standard ---------------
# Result --> set CXXVERSION variable.
# Set ppc_support to TRUE before including this file or ppc and ppc64
# will be treated as invalid architectures since they are no longer supported by Apple
# 
function(detect_cxx_standard output_var)
        file(WRITE "${CMAKE_BINARY_DIR}/cxxversion.cpp" "${cxx_detect_code}")

        enable_language(CXX)

        option(USE_CXX11 "Prefer C++11 features over BOOST whenever possible" OFF)

        if(USE_CXX11)
          try_run(
            run_result_unused
            compile_result_unused
            "${CMAKE_BINARY_DIR}"
            "${CMAKE_BINARY_DIR}/cxxversion.cpp"
            COMPILE_DEFINITIONS "-std=c++11" # Try first with C++11 features enabled
            COMPILE_OUTPUT_VARIABLE CXXVERSION
          )
        endif(USE_CXX11)

        # Parse the architecture name from the compiler output
        string(REGEX MATCH "CXXVERSION ([a-zA-Z0-9_]+)" CXXVERSION "${CXXVERSION}")

        if(NOT CXXVERSION)
          if(USE_CXX11)
            message(FATAL_ERROR "C++11 features were requested (USE_CXX11=ON), but trying -std=c++11 failed.")
          endif(USE_CXX11)

          try_run(
              run_result_unused
              compile_result_unused
              "${CMAKE_BINARY_DIR}"
              "${CMAKE_BINARY_DIR}/cxxversion.cpp"
              COMPILE_OUTPUT_VARIABLE CXXVERSION
          )

          # Parse the architecture name from the compiler output
          string(REGEX MATCH "CXXVERSION ([a-zA-Z0-9_]+)" CXXVERSION "${CXXVERSION}")
        else(NOT CXXVERSION)

          # Turn on C++11 features if they can be specified this way
          SET(CMAKE_CXX_FLAGS "-std=c++11 ${CMAKE_CXX_FLAGS}" PARENT_SCOPE)
        endif(NOT CXXVERSION)

        # Get rid of the value marker leaving just the architecture name
        string(REPLACE "CXXVERSION " "" CXXVERSION "${CXXVERSION}")

        set(${output_var} "${CXXVERSION}" CACHE STRING "C++ standard used to compile the process" FORCE)
endfunction()


# ----- Return the version of the current compiler ----
# Usefull only for cmake version < 3.0
# cf http://stackoverflow.com/questions/435708/any-way-in-cmake-to-require-gcc-version-4
# with cmake version >= 2.8.8 we can use CMAKE_C_COMPILER_VERSION and
# CMAKE_CXX_COMPILER_VERSION
function(detect_cxx_compiler_version gcc_compiler_version)
  if(CMAKE_VERSION VERSION_LESS 3.0.0)
    exec_program(
      ${CMAKE_CXX_COMPILER}
      ARGS                    --version
      OUTPUT_VARIABLE _compiler_output)
    string(REGEX REPLACE ".*([0-9]\\.[0-9]\\.[0-9]).*" "\\1"
      compiler_version ${_compiler_output})
    set(${gcc_compiler_version} ${compiler_version} PARENT_SCOPE)
    message(STATUS "C++ compiler version: ${compiler_version} [${}CMAKE_CXX_COMPILER}]")
  else()
    set(${gcc_compiler_version} ${CMAKE_CXX_COMPILER_VERSION} PARENT_SCOPE)
  endif()
endfunction()

