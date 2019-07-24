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
# Find Lapack libraries and headers

# First : complete (lib and headers) search of blas.
if(NOT BLAS_FOUND)
  find_package(BLASFull)
endif()

if(WITH_BLAS MATCHES matlab)
  set(CLAPACK_HEADER "lapack.h")
  set(LAPACK_INCLUDE_SUFFIXES "extern/include")
  check_lapack_libraries(
    LAPACK_LIBRARIES
    ""
    cheev
    ""
    "mwlapack"
    "${BLAS_LIBRARIES}"
    "")

  # TO BE DONE/CHECKED?
elseif(WITH_BLAS MATCHES mkl)
  set(LAPACK_HEADER mkl_lapacke.h)
elseif(WITH_BLAS MATCHES accelerate)
  set(LAPACK_HEADER clapack.h)
elseif(NOT DEFINED WITH_BLAS AND APPLE)
  # The 'default' case on macosx
  set(LAPACK_HEADER clapack.h)
else()
  set(LAPACK_HEADER lapacke.h)
endif()


message("  __ _ __ _  BLAS  INPUTS with blas : ${WITH_BLAS}")
message("BLAS include : ${BLAS_INCLUDE_SUFFIXES}")
message("BLAS root : ${BLAS_ROOT}")

# Call cmake FindLAPACK to find libraries.
find_package(LAPACK)

if(LAPACK_FOUND)
  message(STATUS "Lapack libraries have been found : ${LAPACK_LIBRARIES}. We now turn to headers.")

  unset(LAPACK_FOUND CACHE)
  message(STATUS "Lapack libraries have been found : ${LAPACK_LIBRARIES}. We now turn to headers.")

  # Get an hint for the location of lapack headers : probably BLAS_LIBRARIES/../../include
  if(NOT LAPACK_LIBRARY_DIR)
    # Get first lib in the list ...
    list(GET LAPACK_LIBRARIES 0 LAPACK_LIB)
    # and its parent directory.
    get_filename_component(LAPACK_LIBRARY_DIR ${LAPACK_LIB} DIRECTORY)
  endif()

  find_path(LAPACK_INCLUDE_DIRS
    NAMES ${LAPACK_HEADER}
    HINTS ${LAPACK_LIBRARIES} ${BLAS_ROOT}
    PATH_SUFFIXES ${BLAS_INCLUDE_SUFFIXES}
    )
  
  
  #       set(LAPACK_LIBRARY_DIR ${_bdir} CACHE PATH "Lapack libraries location." FORCE)
  #     endif()
      
  #     # Either LAPACK_DIR has been set by user or we use LAPACK_LIBRARIES location to set it. 
  #     if(NOT LAPACK_DIR)
  #       get_filename_component(_bdir ${LAPACK_LIBRARY_DIR} PATH)
  #       set(LAPACK_DIR ${_bdir} CACHE PATH "Lapack implementation location." FORCE)
  #     endif()
      
  #     set(LAPACK_INC_DIR ${LAPACK_DIR}/include)
  #     ## Use LAPACK_DIR to set path for cmake search. 
  #     ## If pkg-config has been used, maybe it gives some hints about headers location ...
  #     set(CMAKE_INCLUDE_PATH ${LAPACK_INC_DIR} ${INCLUDE_DIR_HINTS})
  #     set(CMAKE_PREFIX_PATH ${LAPACK_DIR})

  #     ## Note Franck : it seems that find_path process does not work as expected on Macosx : it does not respect the search sequence described
  #     # in cmake doc (i.e. CMAKE_PREFIX, then HINTS, PATHS ... ) and always turn to find apple framework cblas if
  #     # NO_DEFAULT_PATH is not set. 
  #     if(APPLE) # First check in HINTS, no default, then global search.
  #       if(LAPACKE_HEADER)
  #         find_path(LAPACK_INCLUDE_DIRS
  #           NAMES ${LAPACKE_HEADER}
  #           HINTS ${LAPACK_DIR} ${LAPACK_INC_DIR}
  #           PATH_SUFFIXES ${LAPACK_INCLUDE_SUFFIXES}
  #           NO_DEFAULT_PATH
  #           )
  #         find_path(LAPACK_INCLUDE_DIRS
  #           NAMES ${LAPACKE_HEADER}
  #           PATH_SUFFIXES ${LAPACK_INCLUDE_SUFFIXES}
  #           )
  #         if(LAPACK_INCLUDE_DIRS)
  #           set(HAS_LAPACKE 1 CACHE BOOL "Lapacke interface is available.")
  #         endif()

  #       endif() # if lapacke.h not needed or not found search for clapack.h ...
  #       if(NOT LAPACK_INCLUDE_DIRS)
  #         find_path(LAPACK_INCLUDE_DIRS
  #           NAMES ${CLAPACK_HEADER}
  #           HINTS ${LAPACK_LIBRARIES}
  #           PATH_SUFFIXES ${LAPACK_INCLUDE_SUFFIXES}
  #           NO_DEFAULT_PATH
  #           )
  #         find_path(LAPACK_INCLUDE_DIRS
  #           NAMES ${CLAPACK_HEADER}
  #           PATH_SUFFIXES ${LAPACK_INCLUDE_SUFFIXES}
  #           )
  #         if(LAPACK_INCLUDE_DIRS)
  #           set(HAS_CLAPACK 1 CACHE BOOL "clapack interface is available.")
  #         else()
  #           if(LAPACK_FIND_REQUIRED)
  #             message(FATAL_ERROR "Lapack headers not found.")
  #           endif()
  #         endif()
  #       endif()

  #     else(APPLE) # The case which is supposed to always work
  #       if(LAPACKE_HEADER)
  #         find_path(LAPACK_INCLUDE_DIRS 
  #           NAMES ${LAPACKE_HEADER}
  #           HINTS ${LAPACK_DIR} ${LAPACK_INC_DIR}
  #           PATH_SUFFIXES ${LAPACK_INCLUDE_SUFFIXES}
  #           )
  #         if(LAPACK_INCLUDE_DIRS)
  #           set(HAS_LAPACKE 1 CACHE BOOL "Lapacke interface is available.")
  #         endif()
  #       endif()
  #       # if lapacke.h not needed or not found search for clapack.h ...
  #       if(NOT LAPACK_INCLUDE_DIRS)
  #         find_path(LAPACK_INCLUDE_DIRS 
  #           NAMES ${CLAPACK_HEADER}
  #           HINTS ${LAPACK_DIR} ${LAPACK_INC_DIR}
  #           PATH_SUFFIXES ${LAPACK_INCLUDE_SUFFIXES}
  #           )
  #         if(LAPACK_INCLUDE_DIRS)
  #           set(HAS_CLAPACK 1 CACHE BOOL "clapack interface is available.")
  #         endif()

  #       endif()
  #     endif()
  #     if(LAPACK_INCLUDE_DIRS)
  #       set(LAPACK_FOUND 1 CACHE BOOL "Found lapack headers dir.")
  #     endif()
      
  #   endif()
  # endif()

  #  if(LAPACK_FOUND) # SiconosConfig.h setup
  #    if(NOT LAPACK_LIBRARIES MATCHES "lapacke.*") # atlas can be the BLAS provider
  #     ## If lapack was found as "generic" but is part of atlas. 
  #     if(LAPACK_LIBRARIES MATCHES "atlas.*") 
  #        set(WITH_LAPACK "atlas" CACHE STRING "lapack implementation type [mkl/openblas/atlas/accelerate/generic]" FORCE)
  #      endif()
  #    endif()

  #   if(WITH_LAPACK STREQUAL "mkl")
  #     set(HAS_MKL_LAPACKE 1 CACHE BOOL "Blas comes from Intel MKL.")
  #     set(LAPACK_SUFFIX)
  #     set(LAPACK_PREFIX "LAPACKE_")

  #   elseif(WITH_LAPACK STREQUAL "accelerate")
  #     set(HAS_ACCELERATE 1 CACHE BOOL "Blas/Lapack come from Accelerate framework ")
  #     set(LAPACK_SUFFIX "_")
  #     set(LAPACK_PREFIX)

  #   elseif(WITH_LAPACK STREQUAL "MATLAB")
  #     set(HAS_MATLAB_LAPACK 1 CACHE BOOL "Blas/Lapack come from MATLAB ")
  #     set(LAPACK_SUFFIX "_")
  #     set(LAPACK_PREFIX)
  #     set(HAS_LAPACKE FALSE)
  #     set(HAS_CLAPACK FALSE)

  #   elseif(WITH_LAPACK STREQUAL "atlas")
  #     set(HAS_ATLAS_LAPACK 1 CACHE BOOL "Blas  comes from Atlas framework ")
  #     set(LAPACK_SUFFIX)
  #     set(LAPACK_PREFIX "clapack_")
  #     set(HAS_LAPACKE FALSE)

  #   elseif(WITH_LAPACK STREQUAL "openblas")
  #     set(HAS_OpenBLAS_LAPACK 1 CACHE BOOL "Blas/Lapack come from OpenBlas ")
  #     set(LAPACK_SUFFIX)
  #     if(HAS_LAPACKE)
  #       set(LAPACK_PREFIX "LAPACKE_")
  #     else()
  #     endif()
  #   else()
  #     set(HAS_GenericCLAPACK 1 CACHE BOOL "Lapack is available from an unknown version.")
  #     if(HAS_LAPACKE)
  #       set(LAPACK_PREFIX "LAPACKE_")
  #     else()
  #       set(LAPACK_PREFIX "clapack_")
  #     endif()
  #   endif()
  # endif()

  # if(HAS_LAPACKE)
  #   set(LAPACK_HEADER ${LAPACKE_HEADER} CACHE STRING "Lapack header name.")
  # else()
  #   set(LAPACK_HEADER ${CLAPACK_HEADER} CACHE STRING "Lapack header name.")
  # endif()    

  # === Extra Checks for lapack functions ===
  # We check only the functions that are known to be un-implemented
  # in some lapack libraries (namely atlas ...)
  # This is probably a temporary check since it's likely
  # we will stop atlas checking for lapack? 
  ## Check if FuncName is available in lapack lib (i.e in one of LAPACK_LIBRARIES)
  ## and if FuncName symbol is present in file Header. 
  # If both are true, result is true.
  # function(check_clapack_has_function genericName FuncName Header result)

  #   check_function_exists(${FuncName} TEST_FUNC)
  #   check_symbol_exists(${FuncName} ${Header} TEST_HEAD)

  #   if(TEST_HEAD AND TEST_FUNC)
  #     set(${result} 1 CACHE BOOL "${genericName} is available in lapack.")
  #   endif()

  #   unset(TEST_HEAD)
  #   unset(TEST_FUNC)
  # endfunction()


  # set(CMAKE_REQUIRED_LIBRARIES ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
  # set(CMAKE_REQUIRED_INCLUDES ${BLAS_INCLUDE_DIRS} ${LAPACK_INCLUDE_DIRS})
  # if (BLA_STATIC AND CMAKE_COMPILER_IS_GNUG77 AND NOT MSVC)
  #   if (NOT GFORTRAN_LIB)
  #     set(GFORTRAN_LIB "gfortran")
  #   endif()
  #   set(CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES} ${GFORTRAN_LIB})
  # endif()

  # if(MSVC AND HAS_LAPACKE)
  #   SET(CMAKE_REQUIRED_FLAGS ${CMAKE_REQUIRED_FLAGS} "-DHAVE_LAPACK_CONFIG_H -DLAPACK_COMPLEX_STRUCTURE")
  #   SET(CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES} ${EXTRA_LAPACK_LIB})
  # endif(MSVC AND HAS_LAPACKE)

  # if(HAS_MATLAB_LAPACK)
  #   SET(CMAKE_REQUIRED_FLAGS ${CMAKE_REQUIRED_FLAGS} "-Dptrdiff_t=long")
  # endif(HAS_MATLAB_LAPACK)

  # ## dgesvd ##
  # unset(HAS_LAPACK_DGESVD)
  # set(GENERIC_NAME "DGESVD")
  # set(FUNC_NAME "${LAPACK_PREFIX}dgesvd${LAPACK_SUFFIX}")
  # check_clapack_has_function(${GENERIC_NAME} ${FUNC_NAME} ${LAPACK_HEADER} HAS_LAPACK_DGESVD)

  # ## dgels ##
  # unset(HAS_LAPACK_DGELS)
  # set(GENERIC_NAME "DGELS")
  # set(FUNC_NAME "${LAPACK_PREFIX}dgels${LAPACK_SUFFIX}")
  # check_clapack_has_function(${GENERIC_NAME} ${FUNC_NAME} ${LAPACK_HEADER} HAS_LAPACK_DGELS)

  # ## dtrtrs ##
  # unset(HAS_LAPACK_DTRTRS)
  # set(GENERIC_NAME "DTRTRS")
  # set(FUNC_NAME "${LAPACK_PREFIX}dtrtrs${LAPACK_SUFFIX}")
  # check_clapack_has_function(${GENERIC_NAME} ${FUNC_NAME} ${LAPACK_HEADER} HAS_LAPACK_DTRTRS)

  # ## Final settings ...
  # set(LAPACK_LIBRARIES ${LAPACK_LIBRARIES} CACHE INTERNAL "Lapack libraries.")

  # if(NOT LAPACK_FOUND AND LAPACK_FIND_REQUIRED)
  #   message(FATAL_ERROR "Cannot find a library with LAPACK API. Please specify library location using LAPACK_DIR option or set your environment variables properly.")
  # endif (NOT LAPACK_FOUND AND LAPACK_FIND_REQUIRED)
  # if(NOT LAPACK_FIND_QUIETLY)
  #   if(LAPACK_FOUND)
  #     message(STATUS "Found a library with LAPACK API (${WITH_LAPACK}).")
  #     message("    Lapack libraries : ${LAPACK_LIBRARIES}.")
  #     message("    Lapack headers : ${LAPACK_HEADER} in ${LAPACK_INCLUDE_DIRS}.")
  #   else(LAPACK_FOUND)
  #     message(STATUS "Cannot find a library with LAPACK API. Maybe you can try again using LAPACK_DIR option or set your environment variables properly.")
  #   endif(LAPACK_FOUND)
  # endif(NOT LAPACK_FIND_QUIETLY)
  message(" LAPACK HEADER are in ;.. ${LAPACK_INCLUDE_DIRS}")

endif()

# -- Library setup --
find_package_handle_standard_args(LAPACK LAPACK_INCLUDE_DIRS)
if(NOT TARGET LAPACK::LAPACK)
  add_library(LAPACK::LAPACK IMPORTED INTERFACE)
  set_property(TARGET LAPACK::LAPACK PROPERTY INTERFACE_LINK_LIBRARIES ${LAPACK_LIBRARIES})
  if(LAPACK_INCLUDE_DIRS)
    set_target_properties(LAPACK::LAPACK PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES "${LAPACK_INCLUDE_DIRS}")
  endif()
endif()

if(LAPACK_FOUND) # SiconosConfig.h setup
  if(NOT LAPACK_LIBRARIES MATCHES "lapacke.*") # atlas can be the BLAS provider
    ## If lapack was found as "generic" but is part of atlas. 
    if(LAPACK_LIBRARIES MATCHES "atlas.*") 
      set(WITH_LAPACK "atlas" CACHE STRING "lapack implementation type [mkl/openblas/atlas/accelerate/generic]" FORCE)
    endif()
  endif()
  
  if(WITH_LAPACK STREQUAL "mkl")
    set(HAS_MKL_LAPACKE 1 CACHE BOOL "Blas comes from Intel MKL.")
    # set(LAPACK_SUFFIX)
    # set(LAPACK_PREFIX "LAPACKE_")
    
  elseif(WITH_LAPACK STREQUAL "accelerate")
    set(HAS_ACCELERATE 1 CACHE BOOL "Blas/Lapack come from Accelerate framework ")
    # set(LAPACK_SUFFIX "_")
    # set(LAPACK_PREFIX)
    
  elseif(WITH_LAPACK STREQUAL "MATLAB")
    set(HAS_MATLAB_LAPACK 1 CACHE BOOL "Blas/Lapack come from MATLAB ")
    # set(LAPACK_SUFFIX "_")
    # set(LAPACK_PREFIX)
    # set(HAS_LAPACKE FALSE)
    # set(HAS_CLAPACK FALSE)
    
  elseif(WITH_LAPACK STREQUAL "atlas")
    set(HAS_ATLAS_LAPACK 1 CACHE BOOL "Blas  comes from Atlas framework ")
    # set(LAPACK_SUFFIX)
    # set(LAPACK_PREFIX "clapack_")
    # set(HAS_LAPACKE FALSE)
    
  elseif(WITH_LAPACK STREQUAL "openblas")
    # set(HAS_OpenBLAS_LAPACK 1 CACHE BOOL "Blas/Lapack come from OpenBlas ")
    # set(LAPACK_SUFFIX)
    # if(HAS_LAPACKE)
    #   set(LAPACK_PREFIX "LAPACKE_")
    # else()
    # endif()
  else()
    # set(HAS_GenericCLAPACK 1 CACHE BOOL "Lapack is available from an unknown version.")
    # if(HAS_LAPACKE)
    #   set(LAPACK_PREFIX "LAPACKE_")
    # else()
    #   set(LAPACK_PREFIX "clapack_")
    # endif()
  endif()
endif()




