#  Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.
#
# Copyright 2018 INRIA.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# --
# 
# Find BLAS libraries and headers
#
# This module is based on FindBlas provided by cmake
# (check details on https://cmake.org/cmake/help/v3.13/module/FindBLAS.html?highlight=findblas)
# 
# but extend it to find also headers.
#
# Usage :
# 
# find_package(BLASFull REQUIRED)
# target_link_libraries(yourlib PRIVATE BLAS::BLAS)
#
# It will handles both linking and include directories for your target.
#
# This module sets the following variables:
#
# BLAS_FOUND - set to true if a library implementing the BLAS interface
#    is found
#  BLAS_LIBRARIES - uncached list of libraries (using full path name) to
#    link against to use BLAS
#  BLAS_INCLUDE_DIRS - location of blas headers found by cmake
#  BLAS_LIBRARY_DIR - location of the first blas lib found by cmake.
#  BLAS_HEADER - name(s) of the header(s) required for cblas.
#
# To find a specific blas set :
# 
#  WITH_BLAS  if set checks only a specific implementation of blas which may be : mkl, openblas, atlas, accelerate, veclib, generic.
# Or :
#  BLAS_DIR if set, cmake uses this path (absolute) to find a blas implementation in BLAS_DIR/include and/or BLAS_DIR/lib
#
# If neither WITH_BLAS nor BLAS_DIR is set, all implementation are searched in this order : mkl, openblas, accelerate, veclib, atlas, generic.
# If available, pkg-config system is used to give hints to cmake search process.
#

# Three scenarii:
# - set WITH_BLAS to choose a specific vendor
# - set BLAS_DIR to path where your distribution is available
# - Neither BLAS_DIR nor WITH_BLAS : let cmake find something on your os, in standard paths.
# 
#
# - WITH_BLAS=<mkl, openblas, atlas, accelerate, matlab>
#
# if :
# - mkl : hint from MKLROOT
# - openblas :
# - accelerate
# - atlas
# - matlab : hint from matlabroot

# Remark : if BLA_VENDOR is set, WITH_BLAS is ignored and overwritten.

# Check BLA_VENDOR
if(BLA_VENDOR)
  if(BLA_VENDOR MATCHES Intel)
    set(WITH_BLAS mkl)
  elseif(BLA_VENDOR MATCHES Apple)
    set(WITH_BLAS accelerate)
  elseif(BLA_VENDOR MATCHES OpenBLAS)
    set(WITH_BLAS openblas)
  endif()
endif()

# Process user option WITH_BLAS (if any)
if(WITH_BLAS)
  list(LENGTH WITH_BLAS blas_opt_len)
  if(blas_opt_len GREATER 1)
    list(GET WITH_BLAS 1 BLAS_ROOT)
    message("input opt : ${WITH_BLAS} of length ${blas_opt_len} in ${BLAS_ROOT}")
  endif()
endif()

if(WITH_BLAS MATCHES matlab)
  # Matlab case is a bit different
  # Header is blas.h in <matlabroot>/extern/include
  # Libs are in <matlabroot>/bin/...
  # So we need matlabroot and lib path (e.g. with env LIB/DYLD_LIBRARY_PATH/LD_LIBRARY_PATH)
  #
  set(BLAS_HEADER blas.h CACHE STRING "Blas header name")
  set(BLAS_INCLUDE_SUFFIXES "extern/include")
  if(${blas_opt_len} LESS 2)
    message(FATAL_ERROR "\nPlease run cmake with -DWITH_BLAS=\"matlab;<matlabroot>\" where matlabroot is the \
path to matlab install (output from matlabroot command in Matlab).\n")
  endif()
  list(GET WITH_BLAS 1 MATLAB_ROOT)
  message(" Try to find Blas in matlab distribution (${MATLAB_ROOT})")
  
elseif(WITH_BLAS MATCHES mkl)
  if(${blas_opt_len} LESS 2)
    message(FATAL_ERROR "\nPlease run cmake with -DWITH_BLAS=\"mkl;<mklroot>\" where mklroot is the \
path to mkl install.\n")
  endif()
endif()

if(NOT WITH_BLAS MATCHES matlab)
  # Call cmake FindBLAS to find libraries.
  find_package(BLAS)
endif()

if(BLAS_FOUND)
  message(STATUS "Blas libraries have been found : ${BLAS_LIBRARIES}. We now turn to headers.")
endif()

## Now the headers ...
# Steps :
# - set BLAS_HEADER name, this depends on each blas implementation.
# - use BLAS_LIBRARIES to provide hints for possible location of header files.
# - search for BLAS_HEADER
#

# Set BLAS_HEADER, name of the header for cblas interface.
# This is the file that will be searched for
# and its path will be appended to BLAS_INCLUDE_DIRS.
# Default (most common?) blas header is cblas.h
if(WITH_BLAS MATCHES matlab)
  set(BLAS_HEADER blas.h CACHE STRING "Blas header name")
  set(BLAS_INCLUDE_SUFFIXES "extern/include")
elseif(WITH_BLAS MATCHES mkl)
  set(BLAS_HEADER mkl_cblas.h CACHE STRING "Blas header name")
else()
  set(BLAS_HEADER cblas.h CACHE STRING "Blas header name")
  if(WITH_BLAS MATCHES atlas)
    set(BLAS_INCLUDE_SUFFIXES atlas)
  elseif(WITH_BLAS MATCHES openblas)
    set(BLAS_INCLUDE_SUFFIXES "openblas")
  elseif(WITH_BLAS MATCHES accelerate)
    set(BLAS_INCLUDE_SUFFIXES Headers Frameworks)
  elseif(NOT DEFINED WITH_BLAS AND APPLE)
    # The 'default' case on macosx
    set(BLAS_INCLUDE_SUFFIXES Headers Frameworks)
  endif()
endif()

if(BLAS_FOUND) # libraries has been properly found.
  unset(BLAS_FOUND CACHE)
  
  # Get an hint for the location of blas headers : probably BLAS_LIBRARIES/../../include
  # Get first lib in the list ...
  list(GET BLAS_LIBRARIES 0 BLAS_LIB)
  # and its parent directory.
  get_filename_component(BLAS_LIBRARY_DIR ${BLAS_LIB} DIRECTORY)

  # Either BLAS_ROOT has been set (in WITH_BLAS list) or we use BLAS_LIBRARIES location to set it.
  # Most of the time, blas org is like : BLAS_ROOT/lib/XX/libcblas.so.
  if(NOT BLAS_ROOT)
    get_filename_component(_bdir ${BLAS_LIBRARY_DIR} DIRECTORY)
    get_filename_component(_bdir ${_bdir} DIRECTORY)
    get_filename_component(BLAS_ROOT ${_bdir} DIRECTORY)
  endif()


  if(WITH_BLAS STREQUAL "Apple" OR WITH_BLAS STREQUAL "accelerate")
    message("Apple accelerate ...")
    # Apple accelerate framework
    find_path(BLAS_INCLUDE_DIRS
      NAMES ${BLAS_HEADER}
      HINTS ${BLAS_LIBRARIES} ${BLAS_ROOT}
      PATH_SUFFIXES ${BLAS_INCLUDE_SUFFIXES}
      NO_DEFAULT_PATH
      )
  else()
    find_path(BLAS_INCLUDE_DIRS
      NAMES ${BLAS_HEADER}
      HINTS ${BLAS_LIBRARIES} ${BLAS_ROOT}
      PATH_SUFFIXES include ${BLAS_INCLUDE_SUFFIXES}
      )
  endif()

endif()

# -- Library setup --
find_package_handle_standard_args(BLAS BLAS_INCLUDE_DIRS)
if(NOT TARGET BLAS::BLAS)
  add_library(BLAS::BLAS IMPORTED INTERFACE)
  set_property(TARGET BLAS::BLAS PROPERTY INTERFACE_LINK_LIBRARIES ${BLAS_LIBRARIES})
  if(BLAS_INCLUDE_DIRS)
    set_target_properties(BLAS::BLAS PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES "${BLAS_INCLUDE_DIRS}")
  endif()
endif()

# -- Extra setup required for SiconosConfig.h generation. --
if(BLAS_FOUND)
  set(HAS_CBLAS 1 CACHE BOOL "A CBlas implementation is available.")
  
  if(WITH_BLAS MATCHES mkl)
    set(HAS_MKL_CBLAS 1 CACHE BOOL "Blas comes from Intel MKL.")
    
  elseif(WITH_BLAS MATCHES accelerate)
    set(HAS_ACCELERATE 1 CACHE BOOL "Blas/Lapack come from Accelerate framework ")
    
  elseif(WITH_BLAS MATCHES atlas)
    set(HAS_ATLAS_CBLAS 1 CACHE BOOL "Blas  comes from Atlas framework ")
    
  elseif(WITH_BLAS MATCHES openblas)
    set(HAS_OpenBLAS 1 CACHE BOOL "Blas/Lapack come from OpenBlas ")

  elseif(APPLE AND BLAS_LIBRARIES MATCHES Accelerate) # default behavior on macos ...
    set(HAS_ACCELERATE 1 CACHE BOOL "Blas/Lapack come from Accelerate framework ")
  else()
    set(HAS_GenericCBLAS 1 CACHE BOOL "Blas is available from an unknown version.")
  endif()
endif()

