# SPDX-FileCopyrightText: 2026 INRIA
#
# SPDX-License-Identifier: Apache-2.0
#
# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.

#[=======================================================================[.rst:

BlasLapackUtils.cmake :

a set of functions used in Siconos to detect Blas/Lapack implementation
and configure Siconos according to the version found.

#]=======================================================================]

include(CheckFunctionExists)
include(CheckSymbolExists)



# --------------------------------------------
# Try to guess the name of the 'vendor' of the Blas
# implementation that has been found.
# It is based on the list of libraries saved in the variable BLAS_LIBRARIES.
#
# Usage:
# find_package(BLAS)
# guess_blas_name_and_headers()
#
# The following variables will be set in parent scope:
#
# * BLAS_NAME: name of the found implementation
#   (could be: mkl, Matlab, OpenBLAS, Accelerate, Unknown).
# * BLAS_HEADER: name of the required blas header file
# * LAPACK_HEADER: name of the required lapack header file
# * BLAS_PATH_SUFFIXES: dir name(s) used as hints to search for lib or headers
# * LAPACK_EXTRA_LIBS: list of libraries names that may contain the required functions
# * LAPACK_PREFIX, LAPACK_SUFFIX: prefix and suffix for lapack functions names
#   (e.g. in openblas+lapacke : prefix = LAPACKE_)
function(guess_blas_name_and_headers)

  if(BLA_VENDOR MATCHES Intel OR BLAS_LIBRARIES MATCHES mkl)
    # Intel MKL
    set(BLAS_NAME "mkl")
    set(BLAS_HEADER mkl_cblas.h)
    set(LAPACK_HEADER mkl_lapacke.h)
    set(LAPACK_SUFFIX)
    set(LAPACK_PREFIX "LAPACKE_")
  elseif(BLAS_LIBRARIES MATCHES mwblas OR BLA_VENDOR MATCHES "[Mm]atlab")
    # Matlab
    set(BLAS_NAME "Matlab")
    set(BLAS_HEADER blas.h)
    set(LAPACK_HEADER lapack.h)
    set(BLAS_INCLUDE_SUFFIXES "extern/include")
    set(LAPACK_EXTRA_LIBS mwlapack)
    set(LAPACK_SUFFIX "_")
    set(LAPACK_PREFIX)
  elseif(BLAS_LIBRARIES MATCHES openblas OR BLA_VENDOR MATCHES "^OpenBLAS$")
    set(BLAS_NAME "OpenBlas")
    set(BLAS_HEADER cblas.h)
    set(LAPACK_HEADER lapacke.h)
    set(BLAS_INCLUDE_SUFFIXES openblas)
    set(LAPACK_EXTRA_LIBS lapacke)
    set(LAPACK_SUFFIX)
    set(LAPACK_PREFIX "LAPACKE_")
  elseif(BLAS_LIBRARIES MATCHES "Accelerate" OR BLA_VENDOR MATCHES "Apple")
    set(BLAS_NAME "Accelerate")
    set(BLAS_HEADER cblas.h)
    set(LAPACK_HEADER clapack.h)
    set(BLAS_INCLUDE_SUFFIXES Headers Frameworks)
    set(LAPACK_SUFFIX "_")
    set(LAPACK_PREFIX  "")
  else()
    set(BLAS_HEADER cblas.h)
    set(LAPACK_HEADER lapacke.h)
    set(BLAS_NAME "Unknown")
    set(LAPACK_EXTRA_LIBS lapacke)
    set(LAPACK_PREFIX "LAPACKE_")
    set(LAPACK_SUFFIX)
  endif()
  set(BLAS_NAME ${BLAS_NAME} PARENT_SCOPE)
  set(BLAS_HEADER ${BLAS_HEADER} PARENT_SCOPE)
  set(LAPACK_HEADER ${LAPACK_HEADER} PARENT_SCOPE)
  set(BLAS_INCLUDE_SUFFIXES ${BLAS_INCLUDE_SUFFIXES} PARENT_SCOPE)
  set(LAPACK_EXTRA_LIBS ${LAPACK_EXTRA_LIBS} PARENT_SCOPE)
  set(LAPACK_PREFIX ${LAPACK_PREFIX} PARENT_SCOPE)
  set(LAPACK_SUFFIX ${LAPACK_SUFFIX} PARENT_SCOPE)

  message(STATUS "Found a Blas implementation named ${BLAS_NAME}.")
  message(STATUS "Start headers search and functions checking.")

endfunction()

# --------------------------------------------
# Build a list of candidate include directories from a list of
# libraries (<libdir>/../include), so that headers are searched
# next to the library actually found by FindBLAS/FindLAPACK.
#
# Usage: blas_include_hints_from_libraries(OUT_VAR ${BLAS_LIBRARIES})
function(blas_include_hints_from_libraries OUT_VAR)
  set(_hints)
  foreach(_lib IN LISTS ARGN)
    if(IS_ABSOLUTE "${_lib}" AND NOT _lib MATCHES "\\.framework$")
      get_filename_component(_dir "${_lib}" DIRECTORY)
      get_filename_component(_prefix "${_dir}" DIRECTORY)
      list(APPEND _hints "${_prefix}/include" "${_prefix}")
    endif()
  endforeach()
  if(_hints)
    list(REMOVE_DUPLICATES _hints)
  endif()
  set(${OUT_VAR} ${_hints} PARENT_SCOPE)
endfunction()

# --------------------------------------------
# Check if the cblas interface found by cmake
# is complete enough.
# This interface must contain a list of functions (see cblas_functions var below).
#
# Process :
# - check for cblas_dcopy in current BLAS_LIBRARIES.
# - If not found, look for another lib (cblas) and check cblas_dcopy.
# - Finally, in all cases, check for the others required functions
#
# Result : check for a list of function and
# set CBLAS_LIBRARY with some extra libs (if necessary) in parent scope.
function(check_interface NAME)
  set(multiValueArgs FUNCTIONS OPT_FUNCTIONS EXTRA_LIBS)
  set(oneValueArgs PREFIX SUFFIX)
  cmake_parse_arguments(${NAME} "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

  set(CMAKE_REQUIRED_INCLUDES ${${NAME}_INCLUDE_DIR})
  set(CMAKE_REQUIRED_LIBRARIES ${${NAME}_LIBRARIES})
  if(${NAME}DEV_FIND_QUIETLY)
    set(CMAKE_REQUIRED_QUIET TRUE)
  else()
    set(CMAKE_REQUIRED_QUIET FALSE)
  endif()
  # As preconised in CMake doc, we use check_symbol rather than check_function
  # The first function in the list is used :
  list(GET ${NAME}_FUNCTIONS 0 func)
  set(fullname ${${NAME}_PREFIX}${func}${${NAME}_SUFFIX})
  set(CMAKE_REQUIRED_QUIET ON) # No messages
  check_symbol_exists("${fullname}" ${${NAME}_HEADER} HAS_${NAME}_${func})
  # If HAS_${NAME}_${func}, it probably means that the API is OK.
  # Anyway, other required functions will be checked later.
  # Remark: check_symbol_exists save result as an INTERNAL CACHE variable.
  # Else, we check in cblas library :
  if(NOT HAS_${NAME}_${func})
    find_library(C${NAME}_LIBRARY
      NAMES ${${NAME}_EXTRA_LIBS}
      HINTS ${${NAME}_LIBRARIES}
      )
    if(C${NAME}_LIBRARY)
      list(APPEND CMAKE_REQUIRED_LIBRARIES ${C${NAME}_LIBRARY})
      unset(HAS_${NAME}_${func} CACHE)
      check_symbol_exists("${fullname}" ${${NAME}_HEADER} HAS_${NAME}_${func})
      list(APPEND ${NAME}_LIBRARIES ${C${NAME}_LIBRARY})
      set(C${NAME}_LIBRARY ${C${NAME}_LIBRARY} PARENT_SCOPE)
    endif()
    if(NOT HAS_${NAME}_${func})
      if(${NAME}DEV_FIND_REQUIRED)
        message(FATAL_ERROR "${NAME} implementation found is not complete, missing function: ${fullname}.")
      else()
        message(WARNING "${NAME} implementation found is not complete, missing function: ${fullname}.")
      endif()
    endif()
  endif()
  # Extra checks
  set(CMAKE_REQUIRED_LIBRARIES ${${NAME}_LIBRARIES})
  set(CMAKE_REQUIRED_INCLUDES ${${NAME}_INCLUDE_DIR})
  foreach(func IN LISTS ${NAME}_FUNCTIONS)
    set(fullname ${${NAME}_PREFIX}${func}${${NAME}_SUFFIX})
    check_symbol_exists("${fullname}" ${${NAME}_HEADER} HAS_${NAME}_${func})
    if(NOT HAS_${NAME}_${func})
      if(${NAME}DEV_FIND_REQUIRED)
        message(FATAL_ERROR "${NAME} implementation found is not complete, missing function : ${fullname}.")
      else()
        message(WARNING "${NAME} implementation found is not complete, missing function : ${fullname}.")
      endif()
    endif()
  endforeach()
  foreach(func IN LISTS ${NAME}_OPT_FUNCTIONS)
    set(fullname ${${NAME}_PREFIX}${func}${${NAME}_SUFFIX})
    check_symbol_exists("${fullname}" ${${NAME}_HEADER} HAS_${NAME}_${func})
    if(NOT HAS_${NAME}_${func})
      message(STATUS "${NAME} implementation found is not complete, missing function : ${fullname}.")
    endif()
  endforeach()
  set(CMAKE_REQUIRED_LIBRARIES)
  set(CMAKE_REQUIRED_INCLUDES)
  unset(CMAKE_REQUIRED_QUIET)

endfunction()
