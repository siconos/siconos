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


cmake_minimum_required(VERSION 2.8.5)
include(${CMAKE_SOURCE_DIR}/../Numerics/cmake/OutOfSourcesBuild.cmake)

include(CMakeParseArguments)
include(ExternalProject)

# may be changed on the command line
# example : cmake -DINSTALL_COMMAND='gksudo;make;install'
# the ';' are for the separation of each items : command and arguments.

if(NOT INSTALL_COMMAND)
  set(INSTALL_COMMAND sudo make install)
endif()

# default is HEAD
if(NOT GIT_TAG)
  set(GIT_TAG HEAD)
endif()

# default git repo is this cmake source dir
if(NOT GIT_REPOSITORY)
  set(GIT_REPOSITORY ${CMAKE_SOURCE_DIR}/..)
endif()

option(FROM_REPO "Build Siconos from repository. Default = OFF" OFF)
option(WITH_KERNEL "Build Kernel. Default = ON" ON)
option(WITH_MECHANICS "Build Mechanics. Default = off" OFF)
option(WITH_IO "Build IO. Default = off" OFF)
option(WITH_FRONTEND "Build FrontEnd. Default = ON" ON)
option(WITH_EXAMPLES "Build Examples. Default = OFF" OFF) # not tested
option(WITH_DOC "Build Doc. Default = OFF" OFF) # not yet...
option(WITH_TESTING "Build and run tests. Defaul = ON" ON)

if(WITH_TESTING)
  list(APPEND ARGS -DWITH_TESTING=1)
else()
  list(APPEND ARGS -DWITH_TESTING=0)
endif()

if(FROM_REPO)
  # get sources from git
  ExternalProject_Add(Siconos
    GIT_REPOSITORY ${GIT_REPOSITORY}
    GIT_TAG ${GIT_TAG}
    CONFIGURE_COMMAND ${CMAKE_COMMAND} -E echo "no configure"
    BUILD_COMMAND ${CMAKE_COMMAND} -E echo "no build"
    INSTALL_COMMAND ${CMAKE_COMMAND} -E echo "no install"
    )

  set(SOURCE_DIR Siconos-prefix/src/Siconos)
# Note that this force a configure which is done by build if needed :
  set(numerics_depends Siconos) 
else()
  set(SOURCE_DIR ${CMAKE_SOURCE_DIR}/..)
  set(numerics_depends)
endif()

 
macro(AddSiconosProject)
  set(options)
  set(oneValueArgs)
  set(multiValueArgs DEPENDS CMAKE_ARGS)
  cmake_parse_arguments(project "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
  set(_project ${project_UNPARSED_ARGUMENTS})

  ExternalProject_Add(${_project}
    DOWNLOAD_COMMAND ${CMAKE_COMMAND} -E echo "sources are in ${SOURCE_DIR}/${_project}"
    DEPENDS ${project_DEPENDS}
    SOURCE_DIR ${SOURCE_DIR}/${_project}
    CMAKE_ARGS ${project_CMAKE_ARGS}
    INSTALL_COMMAND ${INSTALL_COMMAND}
    TEST_BEFORE_INSTALL ${WITH_TESTING})

  ExternalProject_Add_Step(${_project} forcebuild
    COMMAND ${CMAKE_COMMAND} -E echo "build of ${_project} is forced"
    DEPENDEES configure
    DEPENDERS build
    ALWAYS 1)
endmacro()

AddSiconosProject(Numerics DEPENDS ${numerics_depends} CMAKE_ARGS ${ARGS})

set(frontend_depends Numerics)

if(WITH_KERNEL)
  AddSiconosProject(Kernel DEPENDS Numerics CMAKE_ARGS ${ARGS})
  list(APPEND frontend_depends Kernel)
endif()



if(WITH_MECHANICS)
  AddSiconosProject(Mechanics DEPENDS Kernel CMAKE_ARGS ${ARGS})
  list(APPEND frontend_depends Mechanics)
endif()


if(WITH_IO)
  AddSiconosProject(IO DEPENDS Kernel CMAKE_ARGS ${ARGS})
  list(APPEND frontend_depends IO)
endif()

if(WITH_FRONTEND)
  AddSiconosProject(Front-End DEPENDS ${frontend_depends} CMAKE_ARGS ${ARGS})
endif()

if(WITH_EXAMPLES)
  AddSiconosProject(Examples DEPENDS ${frontend_depends} CMAKE_ARGS ${ARGS})
endif()
