# ============================================================
# Prepare swig config to generate python bindings for siconos
# ============================================================
if(NOT WITH_PYTHON_WRAPPER)
  return()
endif()

find_package(SWIG 3.0 REQUIRED)
include(${SWIG_USE_FILE})

# Name of the generated Python package
set(SICONOS_PYTHON_PACKAGE siconos CACHE INTERNAL "Name of the Siconos python package.")
  
# --------------- Python install setup ---------------
# Set path for siconos-python installation (SICONOS_PYTHON_INSTALL_DIR)
# and get pip install options (PIP_INSTALL_OPTIONS).
include(PythonInstallSetup)
set_python_install_path()
message(STATUS "Siconos python packages will be installed in ${SICONOS_PYTHON_INSTALL_DIR},")
message(STATUS "with the command ${PYTHON_EXECUTABLE} -m pip install ${PIP_INSTALL_OPTIONS} ...")
  
#  ----------- swig -----------

# -- Swig  preprocessor def. common to all swig targets --
# --> SWIG_DEFS

# -- Flags for all swig call ---

# get system architecture 
# https://raw.github.com/petroules/solar-cmake/master/TargetArch.cmake
include(TargetArch)
target_architecture(SYSTEM_ARCHITECTURE)
if(WITH_SYSTEM_INFO) # User defined option, default = off
  include(CMakePrintSystemInformation)
  message(STATUS "SYSTEM ARCHITECTURE: ${SYSTEM_ARCHITECTURE}")
endif()
if(SYSTEM_ARCHITECTURE)
  list(APPEND CMAKE_SWIG_FLAGS "-D__${SYSTEM_ARCHITECTURE}__")
endif()
message(STATUS "SYSTEM ARCHITECTURE: ${CMAKE_SWIG_FLAGS}")

# Generate code with Python 3 specific features and syntax
list(APPEND CMAKE_SWIG_FLAGS "-py3")

# Turn on wrapping of protected members for director classes 
list(APPEND CMAKE_SWIG_FLAGS "-dirprot")
# - without "dirprot", swig wrap public methods and only the 
# protected methods needed to the interface to compile.
# - with "dirprot" swig will attemp to wrap all the public and protected methods at once.

# -dirvtable      - Generate a pseudo virtual table for directors for faster dispatch
# - without "dirprot", swig wrap public methods and only the 
# protected methods needed to the interface to compile.
# - with "dirprot" swig will attemp to wrap all the public and protected methods at once.
set(SWIG_DEFS "-dirprot;-dirvtable;-Wextra")

list(REMOVE_DUPLICATES CMAKE_SWIG_FLAGS)
set(CMAKE_SWIG_FLAGS "${CMAKE_SWIG_FLAGS}" CACHE INTERNAL "Swig flags")

# -- Options related to siconos components
if(WITH_FCLIB)
  set(SWIG_DEFS "${SWIG_DEFS};-DWITH_FCLIB")
endif()
  
if(HAVE_SICONOS_IO)
  set(SWIG_DEFS "${SWIG_DEFS};-DWITH_IO")
  if(WITH_SERIALIZATION)
    set(SWIG_DEFS "${SWIG_DEFS};-DWITH_SERIALIZATION")
  endif()
endif()

# -- Swig files --

# Path to .i files, common to all modules and submodules.
set(SICONOS_SWIG_SOURCE_DIR ${CMAKE_SOURCE_DIR}/wrap/swig
  CACHE INTERNAL "Path to swig files common to all packages.")

# -- Output dir for python generated packages --
set(SICONOS_SWIG_ROOT_DIR ${CMAKE_BINARY_DIR}/wrap/${SICONOS_PYTHON_PACKAGE}
  CACHE INTERNAL "Root path for swig outputs (python packages).")
set(SICONOS_SWIG_BINARY_DIR ${CMAKE_BINARY_DIR}/wrap
  CACHE INTERNAL "Working/binary for swig and python stuff.")

# -- include dir --
set(SICONOS_SWIG_INCLUDE_DIRS ${SICONOS_SWIG_SRC_DIRS} ${SICONOS_SWIG_ROOT_DIR} ${Python3_INCLUDE_DIRS}
  CACHE INTERNAL "Directories required for swig includes.")

foreach(_dir ${SICONOS_SWIG_INCLUDE_DIRS})
  set(SWIG_DEFS "-I${_dir};${SWIG_DEFS}")
endforeach()

set(SWIG_DEFS "${SWIG_DEFS}" CACHE INTERNAL "Swig extra definitions.")

include_directories(${SICONOS_SWIG_INCLUDE_DIRS})

if(WITH_TESTING)
  # windows stuff, probably obsolete ...
  if(CROSSCOMPILING_LINUX_TO_WINDOWS)
    set(EMULATOR "wine")
    set(DRIVE_LETTER "Z:")
  else()
    set(EMULATOR)
    set(DRIVE_LETTER)
  endif()

  # -- test config--
  file(MAKE_DIRECTORY ${SICONOS_SWIG_BINARY_DIR}/tests)
endif()
