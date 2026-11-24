# SPDX-FileCopyrightText: 2026 INRIA
#
# SPDX-License-Identifier: Apache-2.0
#
# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.

#cmake_minimum_required(VERSION 3.2)
#project(
#  ${SKBUILD_PROJECT_NAME}
#  VERSION ${SKBUILD_PROJECT_VERSION}
#  LANGUAGES CXX C)

set(COMPONENT vkernel)
message("-- Set up for ${PROJECT_NAME}_${COMPONENT} library ...\n")
set(${COMPONENT}_DIRS
  src/siconos/vkernel/lib)

# Windows stuff
include(WindowsKernelSetup)

# -- create/setup component target --
include(ComponentSetup)
create_siconos_component(${COMPONENT})

#set(CMAKE_CUDA_FLAGS "-ccbin gcc")
set(CUDA_STANDARD 17)

#set(WITH_GPU 1)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include(CheckCXXCompilerFlag)
set(Siconos_DIR ${SICONOS_INSTALL_DIR})

find_package(Boost REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 3.0.2  CONFIG REQUIRED)

include_directories(${CMAKE_SOURCE_DIR}/externals/boost_hana_pfr)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)

if(WITH_GPU)
  include_directories(${CMAKE_SOURCE_DIR}/externals/cuNSearch/include)
  include_directories(${CMAKE_BINARY_DIR}/externals/cuNSearch/cuNSearch)
else()
  include_directories(${CMAKE_SOURCE_DIR}/externals/CompactNSearch/include)
endif()

# vendored cxd batch narrowphase headers (externals/cxd/)
include_directories(${CMAKE_SOURCE_DIR}/externals/cxd/src)

if(WITH_GPU)
  add_subdirectory(${CMAKE_SOURCE_DIR}/externals/cuNSearch
                   ${CMAKE_BINARY_DIR}/externals/cuNSearch)
else()
  add_subdirectory(${CMAKE_SOURCE_DIR}/externals/CompactNSearch
                   ${CMAKE_BINARY_DIR}/externals/CompactNSearch)
endif()

# every vkernel target includes boost/hana (storage/mp/mp.hpp): outside a
# guix profile (micromamba, brew) the boost headers are not on the default
# include path of the compiler. Not on the vendored neighbour search above:
# in a guix profile the boost include dir is the profile include dir, and
# nvcc given it as -isystem finds the glibc math.h before the libstdc++ one,
# which breaks libcu++ <cuda/std/cmath> (no ::isgreater).
link_libraries(Boost::boost)


link_directories(${CMAKE_BUILD_DIR})

if(WITH_GPU)
  set(flags "-DWITH_GPU")
endif()

set(flags "${flags} -Wall --std=c++23")

# Registry worlds (storage/registry.hpp) check every column lookup, which costs
# about 1% of runtime on bench-disks at 900 disks.  This option drops the
# checks: an access to a column that was never registered becomes undefined
# behaviour instead of an exception.
option(VKERNEL_UNCHECKED_LOOKUP "registry column lookups without checks" OFF)
if(VKERNEL_UNCHECKED_LOOKUP)
  add_compile_definitions(VKERNEL_UNCHECKED_LOOKUP)
endif()
set(fconcepts_diagnostics)
check_cxx_compiler_flag("-fbracket-depth=900" fbracket_depth)
check_cxx_compiler_flag("-fconcepts-diagnostics-depth=1" fconcepts_diagnostics)
check_cxx_compiler_flag("-ftemplate-backtrace-limit=1" ftemplate_backtrace_limit)
check_cxx_compiler_flag("-ftime-trace" ftime_trace)
if(fconcepts_diagnostics)
  set(flags "${flags} -fconcepts-diagnostics-depth=0")
endif()
if(ftemplate_backtrace_limit)
  set(flags "${flags} -ftemplate-backtrace-limit=0")
endif()
if(fbracket_depth)
  set(flags "${flags} -fbracket-depth=900")
endif()

# check_cxx_compiler_flag("-fprofile-instr-generate" fprofile_instr_generate)
# check_cxx_compiler_flag("-fcoverage-mapping" fcoverage_mapping)

# if(fprofile_instr_generate)
#   set(flags "${flags} -fprofile-instr-generate -g")
# endif()
# if(fcoverage_mapping)
#   set(flags "${flags} -fcoverage-mapping")
# endif()

# python_add_library(_nothing MODULE ./src/siconos/pybind/nothing.cpp WITH_SOABI)
# set_target_properties(_nothing PROPERTIES COMPILE_FLAGS "${flags}")
# target_link_libraries(_nothing PRIVATE Eigen3::Eigen Siconos::numerics benchmark pthread fmt pybind11::headers)
# target_compile_definitions(_nothing PRIVATE VERSION_INFO=${PROJECT_VERSION})

# install(TARGETS _nothing DESTINATION nothing)

# The bindings are spread over NONOS_BIND_PARTS translation units (at least
# 2): the simulation alone in the first one, as it instantiates the time step,
# the other items dealt over the others (see part_of in src/siconos/vkernel/py/py_bind.hpp),
# so that they compile in parallel and with less memory each.
set(NONOS_BIND_PARTS 6 CACHE STRING "number of translation units for the _nonos bindings")
set(nonos_bind_sources)
math(EXPR nonos_last_part "${NONOS_BIND_PARTS} - 1")
foreach(part RANGE ${nonos_last_part})
  set(nonos_part_file ${CMAKE_CURRENT_BINARY_DIR}/nonos/bind_part_${part}.cpp)
  file(CONFIGURE OUTPUT ${nonos_part_file}
    CONTENT "#define NONOS_PART ${part}\n#include \"siconos/vkernel/py/py_bind_part.hpp\"\n")
  list(APPEND nonos_bind_sources ${nonos_part_file})
endforeach()

python_add_library(_nonos MODULE
  ./src/siconos/vkernel/py/main.cpp
  ./src/siconos/vkernel/py/world.cpp
  ${nonos_bind_sources}
  WITH_SOABI)

set_target_properties(_nonos PROPERTIES COMPILE_FLAGS "${flags}")
target_compile_definitions(_nonos PRIVATE NONOS_PARTS=${NONOS_BIND_PARTS})

# The binding parts all include py_bind.hpp and recompute the same world-level
# type computations (item lists, attribute transforms, containers): about half
# of their template instantiation time.  In the precompiled header, with
# clang's -fpch-instantiate-templates (added by CMake), they are done once.
# SiconosConfig.h first: the numerics headers depend on its macros (WITH_FCLIB...)
target_precompile_headers(_nonos PRIVATE ${CMAKE_BINARY_DIR}/SiconosConfig.h
                                         ${CMAKE_CURRENT_SOURCE_DIR}/src/siconos/vkernel/py/py_bind.hpp)

# if(ftime_trace)
#   target_compile_options(_nonos PRIVATE
#     -ftime-trace=${CMAKE_CURRENT_BINARY_DIR}/time-trace-nonos.json
#   )
# endif()

if(WITH_GPU)
  target_link_libraries(_nonos PRIVATE mechanics Eigen3::Eigen Boost::boost numerics pthread cuNSearch pybind11::headers)
else()
  target_link_libraries(_nonos PRIVATE mechanics Eigen3::Eigen Boost::boost numerics pthread CompactNSearch pybind11::headers)
endif()


target_compile_definitions(_nonos PRIVATE VERSION_INFO=${PROJECT_VERSION})

set_target_properties(_nonos PROPERTIES
  INSTALL_RPATH "${SICONOS_PYTHON_INSTALL_DIR}/nonos:${CMAKE_INSTALL_PREFIX}/lib64:${CMAKE_INSTALL_PREFIX}/lib"
  BUILD_WITH_INSTALL_RPATH TRUE
  INSTALL_RPATH_USE_LINK_PATH TRUE
)

install(TARGETS _nonos COMPONENT python DESTINATION ${SICONOS_PYTHON_INSTALL_DIR}/nonos/)
install(FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/src/siconos/vkernel/py/nonos/__init__.py
  ${CMAKE_CURRENT_SOURCE_DIR}/src/siconos/vkernel/py/nonos/bridge.py
  ${CMAKE_CURRENT_SOURCE_DIR}/src/siconos/vkernel/py/nonos/simoptions.py
  COMPONENT python
  DESTINATION ${SICONOS_PYTHON_INSTALL_DIR}/nonos/)

add_subdirectory(src)

# # use, i.e. don't skip the full RPATH for the build tree
# set(CMAKE_SKIP_BUILD_RPATH FALSE)

# # when building, don't use the install RPATH already
# # (but later on when installing)
# set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

# set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# # add the automatically determined parts of the RPATH
# # which point to directories outside the build tree to the install RPATH
# set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)


OPTION(BUILD_AS_SHARED_LIBS "Build all the libraries as shared" ON)




# Every public header: the include tree minus py/, tests/, benchmarks/.  The
# trailing slash on src/siconos/ copies the *contents* of that directory, so a
# single statement emits both include/siconos/vkernel/... and
# include/siconos/vkernel.hpp.  Not <component>_INSTALL_INTERFACE_INCLUDE_DIRECTORIES:
# siconos_component_install_setup globs one level and installs flat, which
# would flatten storage/ and simul/ onto one another (both have kernels.hpp).
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/src/siconos/
        DESTINATION include/siconos
        FILES_MATCHING PATTERN "*.hpp"
                       PATTERN "tests" EXCLUDE
                       PATTERN "benchmarks" EXCLUDE
                       PATTERN "py" EXCLUDE)

# The cxd narrowphase (externals/cxd) is vendored but part of the public API:
# collision/cxd_narrowphase.hpp includes its headers, and those include each
# other flat ("cxd_std.hpp") plus a small boost/hana extension they ship
# themselves.  So the directory itself must be on the include path and the
# subtree must be installed with the structure it has in the source tree.
install(DIRECTORY ${CMAKE_SOURCE_DIR}/externals/cxd/src/
        DESTINATION include/siconos/vkernel/cxd
        FILES_MATCHING PATTERN "*.hpp")

# Build root and install root are both "the directory containing siconos/".
# The public headers include Eigen (algebra/eigen_backend.hpp) and several
# Siconos/Numerics headers (CSparseMatrix.h, NumericsMatrix.h, SolverOptions.h,
# SiconosConfig.h, fclib_interface.h), so the component must export both:
# numerics does not re-export Eigen3, and without it the installed headers do
# not compile.  Boost is deliberately not named here: the headers reach it
# through the profile include directory that numerics already exports, while
# siconos-config.cmake does no find_dependency(Boost), so naming Boost::boost
# would break consumers instead of helping them.
target_include_directories(vkernel INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/externals/cxd/src>
  $<INSTALL_INTERFACE:include>
  $<INSTALL_INTERFACE:include/siconos/vkernel/cxd>)
target_link_libraries(vkernel INTERFACE numerics mechanics Eigen3::Eigen)
# the broad phase is part of the public surface: collision/neighborhood.hpp
# includes <CompactNSearch.h> and the symbols live in the library, so a
# consumer using a relation must link it.  cuNSearch (WITH_GPU) is not linked
# here: it is in no export set, and even $<BUILD_INTERFACE:cuNSearch> fails the
# build-tree export(EXPORT siconosTargets); its install rule puts the library
# in the python directory.  In-tree targets link it themselves.
if(NOT WITH_GPU)
  target_link_libraries(vkernel INTERFACE CompactNSearch)
endif()

siconos_component_install_setup(${COMPONENT})
# the public headers use std::ranges::views::zip and <print>: C++23, not 20.
target_compile_features(vkernel PUBLIC cxx_std_23)

include(vkernel_tests.cmake)

# the notes of vkernel/docs/notes, their studies: targets vkernel_study_*
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/docs/notes ${CMAKE_CURRENT_BINARY_DIR}/notes)
