# a meta CMakeLists.txt to build all Siconos components
# in a build directory :
# cmake /your/src/dir/siconos/Build

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

option(FROM_REPO "Build Siconos from repository. Default = OFF" OFF)
option(WITH_KERNEL "Build Kernel. Default = ON" ON)
option(WITH_IO "Build IO. Default = off" OFF)
option(WITH_FRONTEND "Build FrontEnd. Default = ON" ON)
option(WITH_TESTING "Build and run tests. Defaul = ON" ON)

if(WITH_TESTING)
  list(APPEND ARGS -DWITH_TESTING=1)
else()
  list(APPEND ARGS -DWITH_TESTING=0)
endif()

macro(AddSiconosProject)
  set(options)
  set(oneValueArgs)
  set(multiValueArgs DEPENDS CMAKE_ARGS)
  cmake_parse_arguments(project "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
  set(_project ${project_UNPARSED_ARGUMENTS})

  if(FROM_REPO)
# to be completed : need just one repo and then we must set sources inside
#    ExternalProject_Add(${_project}
#      DEPENDS ${project_DEPENDS}
#      GIT_DIRECTORY git+ssh://${REPO_USER}@scm.gforge.inria.fr//gitroot/siconos/siconos.git
#      CMAKE_ARGS ${project_CMAKE_ARGS}
#      INSTALL_COMMAND ${INSTALL_COMMAND}
#      TEST_BEFORE_INSTALL ${WITH_TESTING})
  else()
    ExternalProject_Add(${_project}
      DEPENDS ${project_DEPENDS}
      SOURCE_DIR ${CMAKE_SOURCE_DIR}/../${_project}
      CMAKE_ARGS ${project_CMAKE_ARGS}
      INSTALL_COMMAND ${INSTALL_COMMAND}
      TEST_BEFORE_INSTALL ${WITH_TESTING})
  endif()

  ExternalProject_Add_Step(${_project} forcebuild
    COMMAND ${CMAKE_COMMAND} -E echo "build of ${_project} is forced"
    DEPENDEES configure
    DEPENDERS build
    ALWAYS 1)
endmacro()

AddSiconosProject(Numerics CMAKE_ARGS ${ARGS})

set(frontend_depends Numerics)

if(WITH_KERNEL)
  AddSiconosProject(Kernel DEPENDS Numerics CMAKE_ARGS ${ARGS})
endif()

list(APPEND frontend_depends Kernel)

if(WITH_IO)
  AddSiconosProject(IO DEPENDS Kernel CMAKE_ARGS ${ARGS})
  list(APPEND frontend_depends IO)
endif()

if(WITH_FRONTEND)
  AddSiconosProject(Front-End DEPENDS ${frontend_depends} CMAKE_ARGS ${ARGS})
endif()



