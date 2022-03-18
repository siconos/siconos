# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.
#
# Copyright 2022 INRIA.
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
Check if Bullet is required by Siconos.
If so:
- look for the proper Bullet version (3.17 is the min)
- configure siconos to link with Bullet

The following options control if and how Bullet is used

* WITH_BULLET: activate Bullet and look for it (>=3.17) in the standard paths. Stop if not found.
* Bullet_ROOT=<some_path>: activate Bullet and look for it (>=3.17) in <some_path>. Stop if not found
* BULLET_INSTALL: activate Bullet, download and install the proper version at the same place as Siconos (CMAKE_INSTALL_PREFIX). This
location might be used later as input to BULLET_ROOT.

#]=======================================================================]

function(set_bullet_target)
  if(BULLET_DEFINITIONS)
    string(REPLACE "-D" "" BULLET_DEFINITIONS ${BULLET_DEFINITIONS})
  endif()
  create_target(NAME BULLET::BULLET
    LIBRARIES "${BULLET_LIBRARIES}"
    LIBRARY_DIRS "${BULLET_ROOT_DIR}/${BULLET_LIBRARY_DIRS}"
    INCLUDE_DIRS "${BULLET_INCLUDE_DIRS}"
    COMPILE_DEFINITIONS "${BULLET_DEFINITIONS}")
  # Add bullet headers and libs to the build.
  # first draft ... turn this to private later
  target_link_libraries(${COMPONENT} PUBLIC $<BUILD_INTERFACE:BULLET::BULLET>)
  set(SICONOS_HAS_BULLET TRUE CACHE INTERNAL "True if Bullet API has been found and is activated.")
endfunction()

# Three ways:
  # - iNSTALL_BULLET=ON : use fetchcontent to download and install Bullet as a siconos part;
  # - WITH_BULLET is ON, nothing more: look for bullet, check version and link with siconos-mechanics.
  # - User asks explicitely for a specific (already installed) version of Bullet
  #   by providing Bullet_ROOT on cmake command line.
  #   => find it and check the version

  
# Full config :
# - Download, build and install Bullet
# - Create a target BULLET::BULLET
# - Link mechanics with this target
if(INSTALL_BULLET)
  include(FetchContent)
  message(STATUS "Bullet will be downloaded from github repository and installed as a siconos component")
  # Bullet conf. See what's done in https://github.com/bulletphysics/bullet3/blob/master/build_cmake_pybullet_double.sh
  set(BUILD_PYBULLET ON CACHE INTERNAL "")
  set(BUILD_PYBULLET_NUMPY ON CACHE INTERNAL "")
  if(BULLET_USE_DOUBLE_PRECISION)
    set(USE_DOUBLE_PRECISION ON CACHE INTERNAL "")
  endif()
  set(BT_USE_EGL ON CACHE INTERNAL "")
  set(OpenGL_GL_PREFERENCE GLVND CACHE INTERNAL "")
  # Bullet allows cmake 2, setup a few things to avoid warnings.
  cmake_policy(SET CMP0077 NEW) # option() honors normal variables
  cmake_policy(PUSH)
  cmake_policy(SET CMP0072 NEW) # FindOpenGL prefers GLVND
  cmake_policy(PUSH)
  if(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.20")
    cmake_policy(SET CMP0115 OLD) # 
    cmake_policy(PUSH)
  endif()
  set(FETCHCONTENT_QUIET OFF) # verbose mode for fetchcontent. Comment/uncomment according to your needs.
  FetchContent_Declare(Bullet
    GIT_REPOSITORY    https://github.com/bulletphysics/bullet3.git
    GIT_TAG           ebe1916b90acae8b13cd8c6b637d8327cdc64e94 # hash corresponding to tag 3.1.7
    GIT_SHALLOW TRUE
    UPDATE_DISCONNECTED TRUE # Do not update git repo at each run
    LOG_CONFIGURE TRUE
    LOG_BUILD TRUE
    LOG_INSTALL TRUE
    )
  FetchContent_MakeAvailable(Bullet)
  cmake_policy(POP)
  cmake_policy(POP)
  if(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.20")
    cmake_policy(POP)
  endif()
  include(${bullet_BINARY_DIR}/BulletConfig.cmake)
  set(BULLET_INCLUDE_DIRS $<BUILD_INTERFACE:${bullet_SOURCE_DIR}/src> $<INSTALL_INTERFACE:${bullet_INSTALL_DIR}/${BULLET_INCLUDE_DIRS}>)
  set_bullet_target()
  message(STATUS "Built, installed and used bullet-physics version ${BULLET_VERSION_STRING} in ${BULLET_ROOT_DIR}.")
  if(BULLET_USE_DOUBLE_PRECISION)
    message(STATUS "Bullet use double precision for floating point numbers.")
  endif()

elseif(WITH_BULLET OR Bullet_ROOT)
  # Up to cmake 3.22 there is no way to get Bullet version using find_package standard.
  # It's then mandatory to use the 'config' versionb of find_package.
  # Anyway, find_package is not able to check the version since Bullet does not provide a bullet-config-version or BulettConfigVersion file.
  find_package(Bullet CONFIG REQUIRED)
  include(${Bullet_CONFIG})
 
  if(BULLET_VERSION_STRING VERSION_LESS 3.17)
    set(BULLET_FOUND FALSE)
  endif()

  if(NOT BULLET_FOUND)
    message(FATAL_ERROR "Can not find Bullet in the required version (min 3.17). Please try to install it. \

    Alternately: \

    - run cmake for siconos with Bullet_ROOT equal to the path to Bullet version 3.17 \

    - run cmake for siconos with -DINSTALL_BULLET=ON. Bullet will then be installed in ${CMAKE_INSTALL_PREFIX} and Siconos configured to run with Bullet.")
  endif()

  set(BULLET_INCLUDE_DIRS ${BULLET_ROOT_DIR}/${BULLET_INCLUDE_DIRS})
  set_bullet_target()
  message(STATUS "Found bullet-physics version ${BULLET_VERSION_STRING} in ${BULLET_ROOT_DIR}")
  if(BULLET_USE_DOUBLE_PRECISION)
    message(STATUS "Bullet use double precision for floating point numbers.")
  endif()

endif()

