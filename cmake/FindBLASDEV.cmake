#
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
FindBLASDEV
--------

- Find BLAS library AND headers (Fortran and C)
- Check if the functions required by Siconos are available.

First call find_package(BLAS) from cmake which is limited to fortran library
and then search for extra libraries and headers required for the C/C++ interface.

Input Variables
^^^^^^^^^^^^^^^

``BLA_VENDOR``
  As for the module provided by cmake.
  If set, checks only the specified vendor, if not set checks all the
  possibilities. Defaults to Apple on macOS (Accelerate framework), empty
  (all vendors) elsewhere. Changing it in an existing build directory resets
  all cached Blas/Lapack results. List of vendors valid for Siconos:

  * OpenBLAS
  * Matlab
  * Intel10_32 (intel mkl v10 32 bit)
  * Intel10_64lp (intel mkl v10+ 64 bit, threaded code, lp64 model)
  * Intel10_64lp_seq (intel mkl v10+ 64 bit, sequential code, lp64 model)
  * Intel10_64ilp (intel mkl v10+ 64 bit, threaded code, ilp64 model)
  * Intel10_64ilp_seq (intel mkl v10+ 64 bit, sequential code, ilp64 model)
  * Apple
  * Generic

``BLA_PREFER_PKGCONFIG``
  if set pkg-config will be used to search for a BLAS library first
  and if one is found that is preferred

``BLAS_ROOT``: path to blas install if not standard.

Result Variables
^^^^^^^^^^^^^^^^

Same as https://cmake.org/cmake/help/latest/module/FindBLAS.html
plus a target 'BLAS::BLAS'


Usage :

.. code-block:: cmake

    set(BLA_VENDOR OpenBLAS)
    find_package(BLASDEV)

    target_link_libraries(numerics PRIVATE BLAS::BLAS)


It also defines :

``BLAS_INCLUDE_DIR``
  path to headers required for Blas
``BLAS_HEADER``
  name of the blas/cblas header
``BLAS_NAME``
  name of the vendor/implementation found for blas. Possibilities are :

  * mkl
  * Matlab
  * openblas
  * Accelerate
  * Unknown

#]=======================================================================]

include(BlasLapackUtils)

if(NOT BLAS_ROOT)
  set(BLAS_ROOT $ENV{BLAS_ROOT})
endif()

if(BLAS_ROOT)
  set(_BLAS_SEARCH_OPTS
    "HINTS ${BLAS_ROOT} NO_DEFAULT_PATH")
  set(CMAKE_FIND_ROOT_PATH ${BLAS_ROOT})
endif()

# -- Vendor selection --
# When BLA_VENDOR is not set, cmake FindBLAS tries vendors in a fixed order
# where OpenBLAS comes before Apple: on macOS, any OpenBLAS available in
# the search paths (conda, homebrew ...) silently wins over Accelerate.
# => default to Apple Accelerate on macOS. Use -DBLA_VENDOR=OpenBLAS
# (or All, Intel10_64lp ...) to override.
if(NOT DEFINED BLA_VENDOR AND NOT DEFINED ENV{BLA_VENDOR})
  if(APPLE)
    set(_BLA_VENDOR_DEFAULT Apple)
  else()
    set(_BLA_VENDOR_DEFAULT "")
  endif()
  set(BLA_VENDOR "${_BLA_VENDOR_DEFAULT}" CACHE STRING
    "Blas/Lapack vendor (see cmake FindBLAS doc). E.g. Apple, OpenBLAS, Intel10_64lp, Generic, All. Empty means all.")
  unset(_BLA_VENDOR_DEFAULT)
endif()
if(NOT DEFINED BLA_VENDOR AND DEFINED ENV{BLA_VENDOR})
  set(BLA_VENDOR $ENV{BLA_VENDOR})
endif()

# If the vendor changed since the last configure, forget everything
# cached by previous Blas/Lapack searches (libraries, headers,
# symbol checks, HAS_xxx flags), else results get mixed between vendors.
# (A cached BLAS_INCLUDE_DIR without _SICONOS_BLA_VENDOR_USED means a build
# dir configured before this mechanism existed: reset it as well.)
if(NOT "${_SICONOS_BLA_VENDOR_USED}" STREQUAL "${BLA_VENDOR}"
    AND (DEFINED _SICONOS_BLA_VENDOR_USED OR DEFINED CACHE{BLAS_INCLUDE_DIR}))
  message(STATUS "BLA_VENDOR changed ('${_SICONOS_BLA_VENDOR_USED}' -> '${BLA_VENDOR}'), resetting Blas/Lapack cache entries.")
  get_property(_cache_vars DIRECTORY PROPERTY CACHE_VARIABLES)
  foreach(_var IN LISTS _cache_vars)
    if(_var MATCHES "^C?(BLAS|LAPACK)_" AND NOT _var MATCHES "_ROOT$")
      unset(${_var} CACHE)
    elseif(_var MATCHES "^HAS_(BLAS|LAPACK)_"
        OR _var MATCHES "^HAS_(OpenBLAS|ACCELERATE|MKL_CBLAS|MKL_LAPACKE|MATLAB_LAPACK)$")
      unset(${_var} CACHE)
    endif()
  endforeach()
  unset(_cache_vars)
  unset(_var)
endif()
set(_SICONOS_BLA_VENDOR_USED "${BLA_VENDOR}" CACHE INTERNAL "BLA_VENDOR used for the last Blas/Lapack search.")

# First : search for BLAS and use the find process provided by cmake.
# This will only found libraries, no header, no c-interface.
if(BLASDEV_FIND_REQUIRED)
  find_package(BLAS REQUIRED)
else()
  find_package(BLAS)
endif()

# Try to guess which Blas implementation is used, based on the library name/path.
#  - Set blas header, according to the Blas implementation.
#  - Set all info required later for lapack setup (LAPACK_HEADER, LAPACK_PREFIX ...)
guess_blas_name_and_headers()

# Search blas header
if(BLAS_NAME STREQUAL "Accelerate")
  # Apple accelerate framework
  find_path(BLAS_INCLUDE_DIR
    NAMES ${BLAS_HEADER}
    HINTS ${BLAS_ROOT} ${BLAS_LIBRARIES}
    PATH_SUFFIXES ${BLAS_INCLUDE_SUFFIXES}
    NO_DEFAULT_PATH
    )
elseif(BLAS_NAME STREQUAL "mkl")
  # No need to find headers : we will use -qmkl compilation option
  set(BLAS_INCLUDE_DIR $ENV{MKLROOT}/include)
else()
  # Look for headers next to the library first, and never in macOS
  # frameworks (Accelerate also ships a cblas.h, which must not be
  # mixed with another Blas implementation).
  blas_include_hints_from_libraries(_BLAS_INC_HINTS ${BLAS_LIBRARIES})
  set(_SAVED_CMAKE_FIND_FRAMEWORK ${CMAKE_FIND_FRAMEWORK})
  set(CMAKE_FIND_FRAMEWORK NEVER)
  find_path(BLAS_INCLUDE_DIR
    NAMES ${BLAS_HEADER}
    HINTS ${_BLAS_INC_HINTS}
    ${_BLAS_SEARCH_OPTS}
    PATH_SUFFIXES include ${BLAS_INCLUDE_SUFFIXES}
    )
  set(CMAKE_FIND_FRAMEWORK ${_SAVED_CMAKE_FIND_FRAMEWORK})
  unset(_SAVED_CMAKE_FIND_FRAMEWORK)
endif()
if(NOT BLAS_INCLUDE_DIR)
  if(BLASDEV_FIND_REQUIRED AND NOT BLASDEV_FIND_QUIETLY)
    message(FATAL_ERROR "Can not find blas header, ${BLAS_HEADER}.")
  endif()
endif()
# Check if the blas library found is complete enough.
# Sometimes cblas is already included in blas lib
# (which is the case for libs such as MKL, blas from matlab ...)
set(cblas_functions
  dcopy
  dgemv
  daxpy
  dscal
  dnrm2
  ddot
  dgemm)

check_interface(BLAS
  FUNCTIONS ${cblas_functions}
  EXTRA_LIBS cblas
  PREFIX cblas_)

list(APPEND BLAS_LIBRARIES ${CBLAS_LIBRARY})
list(REMOVE_DUPLICATES BLAS_LIBRARIES)

# -- Library setup --
find_package_handle_standard_args(BLASDEV
   REQUIRED_VARS BLAS_LIBRARIES BLAS_INCLUDE_DIR)
if(NOT TARGET BLAS::BLAS)
  add_library(BLAS::BLAS IMPORTED INTERFACE)
endif()

set_property(TARGET BLAS::BLAS PROPERTY INTERFACE_LINK_LIBRARIES "${BLAS_LIBRARIES}")
if(BLAS_NAME STREQUAL "mkl")
  find_package(MKL CONFIG REQUIRED)
  set_target_properties(BLAS::BLAS PROPERTIES COMPILE_OPTIONS "-qmkl")
  set_property(TARGET BLAS::BLAS PROPERTY INTERFACE_LINK_LIBRARIES MKL::MKL)
endif()

if(BLAS_INCLUDE_DIR AND NOT (BLAS_NAME STREQUAL "mkl"))
  set_target_properties(BLAS::BLAS PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${BLAS_INCLUDE_DIR}")
endif()

get_target_property(CHECK BLAS::BLAS INTERFACE_INCLUDE_DIRECTORIES)

set(CMAKE_FIND_ROOT_PATH)
