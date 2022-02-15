# ============================================================
# Prepare swig config to generate python bindings for siconos
# ============================================================
if(NOT WITH_PYTHON_WRAPPER)
  return()
endif()

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
set(SICONOS_SWIG_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/swig
  CACHE INTERNAL "Path to swig files common to all packages.")

# -- Output dir for python generated packages --
set(SICONOS_SWIG_ROOT_DIR ${CMAKE_CURRENT_BINARY_DIR}/${SICONOS_PYTHON_PACKAGE}
  CACHE INTERNAL "Root path for swig outputs (python packages).")
set(SICONOS_SWIG_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}
  CACHE INTERNAL "Working/binary for swig and python stuff.")

# -- include dir --
set(SICONOS_SWIG_INCLUDE_DIRS ${SICONOS_SWIG_SRC_DIRS} ${SICONOS_SWIG_ROOT_DIR} ${Python3_INCLUDE_DIRS}
  CACHE INTERNAL "Directories required for swig includes.")

foreach(_dir ${SICONOS_SWIG_INCLUDE_DIRS})
  set(SWIG_DEFS "-I${_dir};${SWIG_DEFS}")
endforeach()

set(SWIG_DEFS "${SWIG_DEFS}" CACHE INTERNAL "Swig extra definitions.")

include_directories(${SICONOS_SWIG_INCLUDE_DIRS})

# ------------- Tests -------------
if(WITH_TESTING)
  # windows stuff, probably obsolete ...
  if(CROSSCOMPILING_LINUX_TO_WINDOWS)
    set(EMULATOR "wine")
    set(DRIVE_LETTER "Z:")
  else()
    set(EMULATOR)
    set(DRIVE_LETTER)
  endif()

  # About tox: https://tox.wiki/en/latest/
  # See tox.ini and pyproject.toml files.
  # What for? We can run the test suite with "tox ./wrap in the build dir."
  configure_file(tox.ini.in tox.ini)
  configure_file(conftest.py.in conftest.py) # to set data dir or other useful pytest fixtures.
  configure_file(requirements4tests.txt requirements4tests.txt)

  # -- test config--
  file(MAKE_DIRECTORY ${SICONOS_SWIG_BINARY_DIR}/tests)
endif()

# ====== Generates python files required for build/install process ======

# -- setuptools config --
# doc: https://packaging.python.org/en/latest/tutorials/packaging-projects/
# Required files: setup.cfg, MANIFEST.in, pyproject.toml.

configure_file(setup.cfg.in setup.cfg)
configure_file(pyproject.toml pyproject.toml)
if(EXISTS  ${CMAKE_CURRENT_SOURCE_DIR}/MANIFEST.in)
  configure_file(MANIFEST.in MANIFEST.in)
endif()
# packaging
configure_file(${CMAKE_SOURCE_DIR}/COPYING LICENSE)
configure_file(${CMAKE_SOURCE_DIR}/README.md README.md)


# -- siconos package --
message(STATUS "Generate __init__.py file ...")

configure_file(__init__.py.in ${SICONOS_SWIG_ROOT_DIR}/__init__.py)

# -- siconos.utils submodule --
file(COPY utils  DESTINATION ${SICONOS_SWIG_ROOT_DIR}/)

# ====== Create (and setup) build/install target ======
add_custom_target(python-install
  COMMAND ${PYTHON_EXECUTABLE} -m pip install . ${PIP_INSTALL_OPTIONS}
  VERBATIM USES_TERMINAL
  COMMAND_EXPAND_LISTS
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMENT "build/install siconos package")

# execute python-install when target install is called
install(CODE "execute_process(COMMAND ${CMAKE_MAKE_PROGRAM} python-install WORKING_DIRECTORY \"${CMAKE_CURRENT_BINARY_DIR}\")")

# ====== Create a Target to clean sources (remove .pyc files) and build dir ======
file(GLOB_RECURSE PYCFILES "${CMAKE_SOURCE_DIR}/*.pyc")
add_custom_target(pyclean COMMAND rm -f ${PYCFILES}
  COMMAND make clean
  COMMAND rm -rf ${CMAKE_BINARY_DIR}/build
  COMMENT "Clean sources and build.")

