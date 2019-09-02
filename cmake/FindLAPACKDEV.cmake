# Siconos is a program dedicated to modeling, simulation and control
#  of non smooth dynamical systems.
#  Siconos is a free software; you can redistribute it and/or modify
#  it under the terms of the GNU Lesser General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  Siconos is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Lesser General Public License for more details.
# 
#  You should have received a copy of the GNU Lesser General Public License
#  along with Siconos; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
# 
#  Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr

#[=======================================================================[.rst:
FindLAPACKDEV
--------

- Find LAPACK library AND headers (Fortran and C)
- Check if functions required by Siconos are available.

First call find_package(BLASDEV) the find_package(LAPACK) from cmake which is limited to fortran library
and then search for extra libraries and headers required for the C/C++ interface.

Input Variables
^^^^^^^^^^^^^^^

``BLA_VENDOR``
  As for the module provided by cmake.
  If set, checks only the specified vendor, if not set checks all the
  possibilities.  List of vendors valid for Siconos:

  * OpenBLAS
  * ATLAS
  * Matlab
  * Intel10_32 (intel mkl v10 32 bit)
  * Intel10_64lp (intel mkl v10+ 64 bit, threaded code, lp64 model)
  * Intel10_64lp_seq (intel mkl v10+ 64 bit, sequential code, lp64 model)
  * Intel10_64ilp (intel mkl v10+ 64 bit, threaded code, ilp64 model)
  * Intel10_64ilp_seq (intel mkl v10+ 64 bit, sequential code, ilp64 model)
  * Apple
  * Generic

  ``LAPACK_ROOT``: path to lapack install if not standard. Notice that BLAS_ROOT will be searched
  by default, if set.

Result Variables
^^^^^^^^^^^^^^^^

Same as https://cmake.org/cmake/help/latest/module/FindLAPACK.html
plus a target 'LAPACK::LAPACK' 


Usage :
  
.. code-block:: cmake

    set(BLA_VENDOR OpenBlas)
    # or
    set(BLAS_ROOT /home/username/some_local_install/)
    # or
    set(LAPACK_ROOT /home/username/some_local_install/)

    find_package(LAPACKDEV)

    target_link_libraries(numerics PRIVATE LAPACK::LAPACK)


It also defines :

``LAPACK_INCLUDE_DIR``
  path to headers required for Lapack
``LAPACK_HEADER``
  name of the lapack/clapack header

#]=======================================================================]


message(" ---- Start Blas/Lapack search process ---- ")
if(NOT LAPACK_ROOT)
  set(LAPACK_ROOT $ENV{LAPACK_ROOT})
endif()

if(LAPACKDEV_FIND_REQUIRED)
  set(_LAPACK_SEARCH_OPTS REQUIRED)
endif()

find_package(BLASDEV ${_LAPACK_SEARCH_OPTS})

if(BLAS_ROOT OR LAPACK_ROOT)
  set(_SEARCH_OPTS
    HINTS ${BLAS_ROOT} ${LAPACK_ROOT} NO_DEFAULT_PATH)
  set(CMAKE_FIND_ROOT_PATH ${BLAS_ROOT})
  list(APPEND CMAKE_FIND_ROOT_PATH ${LAPACK_ROOT})
endif()

if(NOT BLAS_NAME STREQUAL "mkl")
  find_package(LAPACK ${_LAPACK_SEARCH_OPTS})
endif()

## Now the headers ...
if(BLASDEV_FOUND AND LAPACK_FOUND)
  if(BLAS_NAME STREQUAL "Accelerate")
    # Apple accelerate framework
    find_path(LAPACK_INCLUDE_DIR
      NAMES ${LAPACK_HEADER}
      HINTS ${LAPACK_LIBRARIES} ${BLAS_ROOT}
      PATH_SUFFIXES ${BLAS_INCLUDE_SUFFIXES}
      NO_DEFAULT_PATH
      )
  else()
    find_path(LAPACK_INCLUDE_DIR
      NAMES ${LAPACK_HEADER}
      ${_SEARCH_OPTS}
      PATH_SUFFIXES include ${BLAS_INCLUDE_SUFFIXES}
      )
  endif()

  # SiconosConfig.h setup
  if(BLAS_NAME STREQUAL "mkl")
    set(HAS_MKL_LAPACKE 1 CACHE BOOL "Blas comes from Intel MKL.")
  elseif(BLAS_NAME STREQUAL "OpenBlas")
    set(HAS_OpenBLAS 1 CACHE BOOL "Blas comes from OpenBLAS.")   
  elseif(BLAS_NAME STREQUAL "Accelerate")
    set(HAS_ACCELERATE 1 CACHE BOOL "Blas/Lapack come from Accelerate framework ")
  elseif(BLAS_NAME STREQUAL "Matlab")
    set(HAS_MATLAB_LAPACK 1 CACHE BOOL "Blas/Lapack come from MATLAB ")
  elseif(BLAS_NAME STREQUAL "ATLAS")
    set(HAS_ATLAS_LAPACK 1 CACHE BOOL "Blas  comes from Atlas framework ")
  endif()

  if(MSVC AND HAS_LAPACKE)
    SET(CMAKE_REQUIRED_FLAGS ${CMAKE_REQUIRED_FLAGS} "-DHAVE_LAPACK_CONFIG_H -DLAPACK_COMPLEX_STRUCTURE")
    SET(CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES} ${EXTRA_LAPACK_LIB})
  endif(MSVC AND HAS_LAPACKE)
  
  if(HAS_MATLAB_LAPACK)
    SET(CMAKE_REQUIRED_FLAGS ${CMAKE_REQUIRED_FLAGS} "-Dptrdiff_t=long")
  endif(HAS_MATLAB_LAPACK)
  
  #check_lapack_interface()
  # Functions that are required
  set(lapack_functions
    dgetrf
    dgetri
    dgesv
    dgetrs
    dpotrf
    )

  # Functions that are optional
  set(opt_functions
    dgesvd
    dgels
    dtrtrs)

  check_interface(LAPACK
    FUNCTIONS ${lapack_functions}
    OPT_FUNCTIONS ${opt_functions}
    EXTRA_LIBS ${LAPACK_EXTRA_LIBS}
    PREFIX ${LAPACK_PREFIX}
    SUFFIX ${LAPACK_SUFFIX})

  
  list(APPEND LAPACK_LIBRARIES ${CLAPACK_LIBRARY})
  list(REMOVE_DUPLICATES LAPACK_LIBRARIES)

endif()

# -- Library setup --
find_package_handle_standard_args(LAPACKDEV
  REQUIRED_VARS LAPACK_LIBRARIES LAPACK_INCLUDE_DIR)

if(NOT TARGET LAPACK::LAPACK)
  add_library(LAPACK::LAPACK IMPORTED INTERFACE)
  set_property(TARGET LAPACK::LAPACK PROPERTY INTERFACE_LINK_LIBRARIES "BLAS::BLAS;${LAPACK_LIBRARIES}")
  if(LAPACK_INCLUDE_DIRS)
    set_target_properties(LAPACK::LAPACK PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES "${LAPACK_INCLUDE_DIRS}")
  endif()
endif()

set(CMAKE_FIND_ROOT_PATH)
message(" ---- End of Blas/Lapack search process ---- ")

