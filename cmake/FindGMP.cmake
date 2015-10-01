# Find the Siconos Numerics includes and libraries.
# The following variables are set if gmp library is found.  If gmp is not
# found, GMP_FOUND is set to false.
#  GMP_FOUND        - True when the GMP include directory is found.
#  GMP_INCLUDE_DIR - the path to where the gmp include files are.
#  GMP_LIBRARY_DIRS - The path to where the gmp library files are.
#  GMP_LIBRARIES    - The libraries to link against gmp

# Set GMP_DIR=where_fftw_is_installed if it's not in a "classic" place or if you want a specific version
#
# Note FP : see http://www.cmake.org/Wiki/CMake:How_To_Find_Libraries
INCLUDE(FindPackageHandleStandardArgs)
find_package(PkgConfig)
include(LibFindMacros)

# Use pkg-config to get hints about paths
pkg_check_modules(GMP_PKGCONF QUIET gmp)

# First search library and header in GMP_DIR, if given
if(GMP_DIR)
  find_library(GMP_LIBRARY 
    NAMES gmp 
    PATHS ${GMP_DIR}
    PATH_SUFFIXES lib
    NO_DEFAULT_PATH
    )
  find_path(GMP_INCLUDE_DIR
    NAMES gmp.h
    PATHS ${GMP_DIR}
    PATH_SUFFIXES include
    NO_DEFAULT_PATH)
# else check env variables 
else()
  find_library(GMP_LIBRARY
    NAMES gmp 
    PATHS ENV DYLD_LIBRARY_PATH ENV LD_LIBRARY_PATH ${GMP_PKGCONF_LIBRARY_DIRS}
    NO_DEFAULT_PATH
    )
  if(GMP_LIBRARY) 
    get_filename_component(GMP_LIBRARY_DIR ${GMP_LIBRARY} PATH)
  endif()
  find_path(GMP_INCLUDE_DIR 
    NAMES gmp.h 
    PATHS ENV INCLUDE ${GMP_LIBRARY_DIR} ${GMP_PKGCONF_INCLUDE_DIRS}
    PATH_SUFFIXES include
    NO_DEFAULT_PATH
    )
endif()

# Last try if nothing has been found before
find_library(GMP_LIBRARY NAMES gmp)
if(GMP_LIBRARY)
  get_filename_component(GMP_LIBRARY_DIR "${GMP_LIBRARY}" PATH)
  get_filename_component(GMP_LIBRARY_DIR_DIR "${GMP_LIBRARY_DIR}" PATH)
endif()
find_path(GMP_INCLUDE_DIR NAMES gmp.h HINTS ${GMP_LIBRARY_DIR_DIR}/include)

# Set the required variables 
set(GMP_PROCESS_INCLUDES GMP_INCLUDE_DIR)
set(GMP_PROCESS_LIBS GMP_LIBRARY)
# This will set properly ${GMP_INCLUDE_DIR} ${GMP_LIBRARY} ${GMP_INCLUDE_DIRS} ${GMP_LIBRARIES}
libfind_process(GMP)

FIND_PACKAGE_HANDLE_STANDARD_ARGS(GMP
  REQUIRED_VARS GMP_LIBRARY)
