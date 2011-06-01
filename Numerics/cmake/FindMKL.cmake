# - Find MKL library (intel implementation for lapack and blas, among other things)
#
# This module sets the following variables:
#  MKL_FOUND - set to true if mkl library was found
#    is found
#  MKL_LINKER_FLAGS - uncached list of required linker flags (excluding -l
#    and -L).
#  MKL_LIBRARIES - uncached list of libraries (using full path name) to 
#    link against to use mkl
#  MKL_INCLUDE_DIRS - list of include directories to find  mkl headers

include(CheckFunctionExists)
include(LibFindMacros)

## Look for mkl_core library and headers
# First look in MKLROOT env var, if set or in user defined var MKLDIR
find_library(MKL_LIBRARY 
  NAMES mkl_core
  PATHS ${MKLDIR} ENV MKLROOT 
  PATH_SUFFIXES lib lib/intel64
  )
find_library(MKL_seq
  NAMES mkl_sequential
  PATHS ${MKLDIR} ENV MKLROOT 
  PATH_SUFFIXES lib lib/intel64
  )

set(MKL_EXTRA_LIBRARIES ${MKL_seq})

find_library(MKL_intel
  NAMES mkl_intel_lp64
  PATHS ${MKLDIR} ENV MKLROOT 
  PATH_SUFFIXES lib lib/intel64
  )

set(MKL_EXTRA_LIBRARIES ${MKL_EXTRA_LIBRARIES} ${MKL_intel})

find_path(
  MKL_INCLUDE_DIR
  NAMES mkl_lapacke.h
  PATHS ${MKLDIR} ENV MKLROOT 
  PATH_SUFFIXES include include/intel64
  )

set(MKL_PROCESS_INCLUDES MKL_INCLUDE_DIR)
set(MKL_PROCESS_LIBS MKL_LIBRARY MKL_EXTRA_LIBRARIES)
libfind_process(MKL)

set(CMAKE_REQUIRED_LIBRARIES ${MKL_LIBRARIES} )
check_function_exists(LAPACKE_dgels LAPACK_mkl_works)
set(CMAKE_REQUIRED_LIBRARIES)
if(LAPACK_mkl_works)
  set(COMPLETE_LAPACK_LIBRARIES ${MKL_LIBRARIES})
endif()


