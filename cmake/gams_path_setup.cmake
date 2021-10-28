#  Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.
#
# Copyright 2021 INRIA.
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
#[=======================================================================[.rst:
What is done in this file ?

* if GAMS_ROOT is set, checks for gams (http://www.gams.com) in GAMS_ROOT directory
   -> list gams source files and add them to numerics build
   -> install gams models files 
    
* Look for Path Ferris solver (optional hint : PathFerris_ROOT),
if found, link numerics with Path Ferris.

* Look for Path VI solver (optional hint : PathVI_ROOT),
if found, link numerics with Path Ferris.


Remarks :

* PathVI and PathFerris are only used (libs and headers) at build time
(through PRIVATE option of target_link_libraries).
Note FP : maybe we should propagate them at runtime ?
(if so path headers from externals must be installed)

* Path headers are distributed in Siconos in externals component.

* If GAMS is setup and found, and if PathFerris/VI_ROOT are not set, GAMS_ROOT is used as a hint to find Path headers.




Example of usage :


cmake -DGAMS_ROOT=<path-to-gams-install> -DPathVI_ROOT=... -DPathFerris_ROOT=...

#]=======================================================================]

# --- GAMS ---

if(GAMS_ROOT)

  # Set the list of source files required to build gams
  set(FILES_TO_CHECK "idxcc.c;optcc.c;gamsxcc.c;gmomcc.c;gclgms.h;gamsxcc.h;idxcc.h;optcc.h;gevmcc.h;gmomcc.h")
  set(GAMS_C_API_DIR "apifiles/C/api")
  foreach(_F IN LISTS FILES_TO_CHECK)
    set(_FF ${GAMS_ROOT}/${GAMS_C_API_DIR}/${_F})
    if(EXISTS ${_FF})
      list(APPEND GAMS_C_API_FILES ${_FF})
    else()
      message(FATAL_ERROR "The GAMS C API file ${_F} has not been found in ${GAMS_ROOT}/${GAMS_C_API_DIR}.")
    endif()
  endforeach()

  set(HAVE_GAMS_C_API TRUE BOOL "The GAMS C API has been found")
  
  include(CheckCCompilerFlag)
  check_c_compiler_flag("-Werror=conversion" C_HAVE_WERR_CONV)
  if(C_HAVE_WERR_CONV)
    set_source_files_properties(${GAMS_C_API_FILES} PROPERTIES COMPILE_FLAGS "-Wno-error=conversion")
  endif()
  
  # Add sources to numerics build
  target_sources(numerics PRIVATE ${GAMS_C_API_FILES})
  
  # Add include only at build time
  # PUBLIC at build time to propagate gams include to test.
  # Should we use them at runtime ? 
  target_include_directories(numerics PUBLIC $<BUILD_INTERFACE:${GAMS_ROOT}/${GAMS_C_API_DIR}>)
  target_include_directories(numerics PRIVATE "${GAMS_ROOT}/testlib_ml") # XXX Hack -- xhub
  
  # path to gams model : used in numerics tests, set in SiconosConfig.h
  # --> path to source hardcoded in SiconosConfig.h, bad isn't it ?
  set(GAMS_MODELS_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/share/gams" CACHE FILEPATH "gams")
  
  # path to gams models to be used at runtime
  set(GAMS_MODELS_SHARE_DIR  "${CMAKE_INSTALL_PREFIX}/share/${PROJECT_NAME}/gams" CACHE FILEPATH "gams")
  
  # Install gams files required at runtime
  install(DIRECTORY ${GAMS_MODELS_SOURCE_DIR} DESTINATION share/${PROJECT_NAME})

  # If gams is found, it may be used as hint to search Path solvers.
  if(NOT PathFerris_ROOT)
    set(PathFerris_ROOT  ${GAMS_ROOT})
  endif()
  
  if(NOT PathFerris_ROOT)
    set(PathFerris_ROOT  ${GAMS_ROOT})
  endif()

  # GAMS distributes CPLEX ...
  # if(NOT CPLEX_ROOT)
  #   set(CPLEX_ROOT ${GAMS_ROOT})
  #   find_package(CPLEX)
  #   target_link_libraries(numerics PRIVATE CPLEX::CPLEX)
  # endif()

endif()

# --- Path Ferris --
find_package(PathFerris)
if(PathFerris_FOUND)
  set(HAVE_PATHFERRIS TRUE CACHE INTERNAL "True if Path solver has been found and is activated.")
  target_link_libraries(numerics PRIVATE PathFerris::PathFerris)
  # Notice that path headers are in externals/PATH_SDK/include !
  target_include_directories(numerics PRIVATE $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/externals>)
  # Should we install these headers ?
  # Should we propagate these includes at runtime ?
endif()

# --- Path VI ---
find_package(PathVI)
if(PathVI_FOUND)
  set(HAVE_PATHVI TRUE CACHE INTERNAL "True if Path solver has been found and is activated.")
  target_link_libraries(numerics PRIVATE PathVI::PathVI)
  # Do we need headers from externals/PATH_SDK/include ?
  # target_include_libraries(numerics PRIVATE $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/externals)
  # Should we install these headers ?
  # Should we propagate these includes at runtime ?

  
  # # XXX hack for now ...
  if(HAVE_GAMS_C_API)
    set(HAVE_GAMS_PATHVI TRUE)
  else()
    set(HAVE_GAMS_PATHVI FALSE)
  endif()
  
endif()
