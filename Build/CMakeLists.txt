# a meta CMakeLists.txt to build all Siconos components
# in a build directory :
# cmake /your/src/dir/siconos/Build

# the build may also be done from a repository on a git tag
# default your source directory : 
# cmake -DFROM_REPO=1 

# directly from gforge repo
# cmake -DFROM_REPO=1 -DGIT_REPOSITORY=git+ssh://<yourlogin>@scm.gforge.inria.fr/gitroot/siconos/siconos.git

# build from a tag, examples : 

# cmake -DFROM_REPO=1 -DGIT_TAG=a5912c7
# cmake -DFROM_REPO=1 -DGIT_TAG=3.5.x

# cmake ~/S/Build  -DFROM_REPO=1 -DGIT_REPOSITORY=/home/maurice/src/git/siconos -DINSTALL_COMMAND='gksudo;make;install' -DON_DASHBOARD=1

# CTest options : -DCTEST_OPTIONS='-j2'

# install prefix : -DINSTALL_PREFIX='/tmp'

# SOURCE_DIRECTORY : -DSOURCE_DIRECTORY=/some/where

# Set minimum version for cmake
cmake_minimum_required(VERSION 2.8.7)

# MACOSX_RPATH required from cmake 3.0.
# See http://www.kitware.com/blog/home/post/510
if(POLICY CMP0042 OLD)
  cmake_policy(SET CMP0042 OLD)
endif()
include(${CMAKE_SOURCE_DIR}/../Numerics/cmake/OutOfSourcesBuild.cmake)

include(CMakeParseArguments)
include(ExternalProject)

# may be changed on the command line
# example : cmake -DINSTALL_COMMAND='gksudo;make;install'
# the ';' are for the separation of each items : command and arguments.

if(NOT INSTALL_COMMAND)
  set(INSTALL_COMMAND sudo ${CMAKE_COMMAND} --build . --target install)
endif()

# default is HEAD
if(NOT GIT_TAG)
  set(GIT_TAG master)
endif()

# default git repo is this cmake source dir
if(NOT GIT_REPOSITORY)
  set(GIT_REPOSITORY ${CMAKE_SOURCE_DIR}/..)
endif()

option(FROM_REPO "Build Siconos from repository. Default = OFF" OFF)
option(WITH_KERNEL "Build Kernel. Default = ON" ON)
option(WITH_CONTROL "Build the Control module. Default = off" OFF)
option(WITH_MECHANICS "Build Mechanics. Default = off" OFF)
option(WITH_IO "Build IO. Default = off" OFF)
option(WITH_FRONTEND "Build FrontEnd. Default = ON" ON)
option(WITH_EXAMPLES "Build Examples. Default = OFF" OFF) # not tested
option(WITH_DOCUMENTATION "Build Documentation. Default = OFF" OFF) # 
option(ONLY_DOCUMENTATION "Build only the documentation. Default = OFF" OFF)
option(WITH_TESTING "Build and run tests. Default = ON" ON)
option(ON_DASHBOARD "Send results on a dasboard. Default = OFF" OFF)
option(WITH_DEBUG "Debug mode. Default = OFF" OFF)

if(WITH_DEBUG)
  set(CTEST_OPTIONS "${CTEST_OPTIONS} --debug")
endif()

if(WITH_TESTING)
  list(APPEND ARGS -DWITH_TESTING=1)
else()
  list(APPEND ARGS -DWITH_TESTING=0)
endif()

if(WITH_DOCUMENTATION)
  list(APPEND ARGS -DWITH_DOCUMENTATION=1)
else()
  list(APPEND ARGS -DWITH_DOCUMENTATION=0)
endif()

if(INSTALL_PREFIX)
  set(_project_CMAKE_ARGS ${project_CMAKE_ARGS} -DCMAKE_INSTALL_PREFIX=${INSTALL_PREFIX})
  list(APPEND ARGS -DCMAKE_INSTALL_PREFIX=${INSTALL_PREFIX})
endif()

if(FROM_REPO)

  if(NOT SOURCE_DIR)
    set(SOURCE_DIR ${CMAKE_BINARY_DIR}/Siconos-prefix/src/Siconos)
  endif()
  # get sources from git
  ExternalProject_Add(Siconos
    GIT_REPOSITORY ${GIT_REPOSITORY}
    GIT_TAG ${GIT_TAG}
    CONFIGURE_COMMAND ${CMAKE_COMMAND} -E echo "no configure"
    BUILD_COMMAND ${CMAKE_COMMAND} -E echo "no build"
    INSTALL_COMMAND ${CMAKE_COMMAND} -E echo "no install"
    SOURCE_DIR ${SOURCE_DIR}
    )

# Note that this force a configure which is done by build if needed :
  set(numerics_depends Siconos) 
else()
  if(NOT SOURCE_DIR)
    set(SOURCE_DIR ${CMAKE_SOURCE_DIR}/..)
  endif()
  set(numerics_depends)
endif()

 
macro(AddSiconosProject)
  set(options)
  set(oneValueArgs TEST_BEFORE_INSTALL TEST_AFTER_INSTALL)
  set(multiValueArgs DEPENDS CMAKE_ARGS)
  cmake_parse_arguments(project "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
  set(_project ${project_UNPARSED_ARGUMENTS})
  
  if(NOT project_TEST_BEFORE_INSTALL)
    if(NOT project_TEST_AFTER_INSTALL)
      set(project_TEST_BEFORE_INSTALL 1)
    else()
      set(project_TEST_BEFORE_INSTALL 0)
    endif()
  endif()
  if(NOT project_TEST_AFTER_INSTALL)
    set(project_TEST_AFTER_INSTALL 0)
  endif()

  # do crazy thing because cmake is stubborn
  set(project_CMAKE_ARGS_sanitized)
  foreach(elt ${project_CMAKE_ARGS})
    set(project_CMAKE_ARGS_sanitized ${elt}__00__${project_CMAKE_ARGS_sanitized})
  endforeach(elt)
  STRING(REGEX REPLACE ":" "__22__" project_CMAKE_ARGS_sanitized "${project_CMAKE_ARGS_sanitized}")

  message(STATUS "Adding ${_project}")
  message(STATUS "Sources are in ${SOURCE_DIR}/${_project}")
  message(STATUS "Binaries are in ${CMAKE_BINARY_DIR}/${_project}")

  SET(THE_GENERATOR "${CMAKE_GENERATOR}")

  # Note ONLY_DOCUMENTATION => without dashboard (need a cmake . -DON_DASHBOARD=0 before)
  if(ON_DASHBOARD)
    ExternalProject_Add(${_project}
      DOWNLOAD_COMMAND ${CMAKE_COMMAND} -E echo "sources are in ${SOURCE_DIR}/${_project}"
      CONFIGURE_COMMAND ${CMAKE_COMMAND} -E echo "configuration is done by ctest"
      BUILD_COMMAND ${CMAKE_CTEST_COMMAND} ${CTEST_OPTIONS}
      -DMODULE=${_project}
      -DMODE=${MODE}
      -DCTEST_CMAKE_GENERATOR=${THE_GENERATOR}
      -DSOURCE_DIRECTORY=${SOURCE_DIR}
      -DBINARY_DIRECTORY=${CMAKE_BINARY_DIR}/${_project}
      -DEXTRA_CMAKE_ARGS=${project_CMAKE_ARGS_sanitized}
      -DFROM_REPO=${FROM_REPO}
      -DBUILD_TYPE=${BUILD_TYPE}
      -DWITH_TESTING=${WITH_TESTING}
      -DCROSSCOMPILING_LINUX_TO_WINDOWS=${CROSSCOMPILING_LINUX_TO_WINDOWS}
      -DBUILDNAME_OPTIONS=${BUILDNAME_OPTIONS}
      -S ${SOURCE_DIR}/Build/CTestDriver.cmake
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
      BINARY_DIR ${CMAKE_BINARY_DIR}/${_project}
      DEPENDS ${project_DEPENDS}
      SOURCE_DIR ${SOURCE_DIR}/${_project}
      INSTALL_DIR ${INSTALL_PREFIX}
      INSTALL_COMMAND ${INSTALL_COMMAND}
      TEST_BEFORE_INSTALL ${project_TEST_BEFORE_INSTALL}
      TEST_AFTER_INSTALL ${project_TEST_AFTER_INSTALL})
  else()
    message(STATUS "configuring ${_project}")
    if (WITH_DOCUMENTATION)
      set(_build_target --target doc) # the effective build step is required by install
    endif()
    if (ONLY_DOCUMENTATION)
      set(_install_command ${CMAKE_COMMAND} -E echo "documentation is installed under GeneratedDocs/")
    else()
      set(_install_command ${INSTALL_COMMAND})
    endif()
    ExternalProject_Add(${_project}
      DOWNLOAD_COMMAND ${CMAKE_COMMAND} -E echo "sources are in ${SOURCE_DIR}/${_project}"
      BUILD_COMMAND  ${CMAKE_COMMAND} --build . ${_build_target}
      DEPENDS ${project_DEPENDS}
      SOURCE_DIR ${SOURCE_DIR}/${_project}
      BINARY_DIR ${CMAKE_BINARY_DIR}/${_project}
      CMAKE_ARGS ${project_CMAKE_ARGS}
      INSTALL_DIR ${INSTALL_PREFIX}
      INSTALL_COMMAND ${_install_command}
      TEST_BEFORE_INSTALL ${WITH_TESTING})
  endif()

  ExternalProject_Add_Step(${_project} forcebuild
    COMMAND ${CMAKE_COMMAND} -E echo "build of ${_project} is forced"
    DEPENDEES configure
    DEPENDERS build
    ALWAYS 1)
  
endmacro()

if(WITH_DOCUMENTATION)
  # Add Docs module
  AddSiconosProject(Docs CMAKE_ARGS ${ARGS})
endif()


AddSiconosProject(Numerics DEPENDS ${numerics_depends} CMAKE_ARGS ${ARGS})

set(Siconos_depends Numerics)

if(WITH_KERNEL)
  AddSiconosProject(Kernel DEPENDS ${Siconos_depends} CMAKE_ARGS ${ARGS})
  list(APPEND Siconos_depends Kernel)
endif()

if(WITH_MECHANICS)
  AddSiconosProject(Mechanics DEPENDS ${Siconos_depends} CMAKE_ARGS ${ARGS})
  list(APPEND Siconos_depends Mechanics)
endif()

if(WITH_CONTROL)
  AddSiconosProject(Control DEPENDS ${Siconos_depends} CMAKE_ARGS ${ARGS})
  list(APPEND Siconos_depends Control)
endif()

if(WITH_IO)
  AddSiconosProject(IO DEPENDS ${Siconos_depends} CMAKE_ARGS ${ARGS})
  list(APPEND Siconos_depends IO)
endif()

if(WITH_FRONTEND)
  AddSiconosProject(Front-End DEPENDS ${Siconos_depends} CMAKE_ARGS ${ARGS})
  list(APPEND Siconos_depends Front-End)
endif()

if(WITH_EXAMPLES)
  AddSiconosProject(Examples 
    DEPENDS ${Siconos_depends} 
    CMAKE_ARGS ${ARGS}
    TEST_AFTER_INSTALL 1)
endif()
