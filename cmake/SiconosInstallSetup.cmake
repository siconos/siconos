#  Siconos is a program dedicated to modeling, simulation and control
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
#[=======================================================================[.rst
set_install_path
------------------
 --- Function to compute the install path/options ---

 Usage :
 set_install_path()

 Sets
 - CMAKE_INSTALL_PREFIX
 - PIP_INSTALL_OPTIONS (cache) (useful only if PYTHON_WRAPPER is ON)
 - SICONOS_PYTHON_INSTALL_DIR (cache) with the full path where 
   Siconos python packages will be installed (useful only if PYTHON_WRAPPER is ON)

    
#]=======================================================================]
  

# -- Set install behavior --
# Two things are to be taken into account:
# - install of libraries, includes ... with make install target
# - install of python packages
# This is quite a mess ...
# Standard user should use default behavior and we hope that they use a virtual env for Python (conda, venv ...)
# 
# 1- We forbid CMAKE_INSTALL_PREFIX set by user
# 
# 2- Default leads to "user" install (no root privileges)
#     - make install and python install in CONDA_PREFIX or VIRTUAL_ENV if they exist.
#     - else, make install in $HOME/.siconos and python as if run with pip install --user.
# 
# 3- if SICONOS_INSTALL_SYSTEM_WIDE is True
#     - install with root privileges
#     - in std path (/usr/local ...)
#
# 4- if SICONOS_CUSTOM_INSTALL=somewhere (advanced usage !)
#     - make install in somewhere
#     - if ISOLATED_INSTALL=False (default) install python packages the same way as in the default case above
#     - if ISOLATED_INSTALL=True, install everything (including python packages) in somewhere (this case is useful for guix ...)

# Finding where pip or python will put things is a real nightmare, especially on Debian
# like systems where pip use site-packages, the distro dist-packages and so on.
# See https://www.python.org/dev/peps/pep-0668/
# Note : witth python 3.10, --prefix option of pip always add "local" to install path.
# 
# Thus, we choose:
#
# - default behavior: install in CONDA_PREFIX or VIRTUAL_ENV if they exist. If not run pip install --user
# - system wide install behavior: run pip install, assuming that the caller is root. Else ... well, a warning has been send by cmake earlier.
# - fully isolated install: run pip install --prefix, assuming that the user knows what he's doing.
# 
# Notice that the case: cmake as user in python venv + system wide install as root (outside venv) will lead to unexpected and unwanted results ...
# 
function(set_install_path)

  # --- Ensure that CMAKE_INSTALL_PREFIX is not set by user
  if(NOT CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    if(NOT EXISTS ${CMAKE_BINARY_DIR}/CMakeCache.txt) 
      # if CMAKE_INSTALL_PREFIX has been explicitely set and if it's the first cmake run
      message(FATAL_ERROR "Please do not set CMAKE_INSTALL_PREFIX.")
    endif()
  endif()

  # Check python env
  set(PIP_INSTALL_OPTIONS_LOCAL)
  execute_process(COMMAND ${Python_EXECUTABLE} -c "import sys;print(sys.prefix != sys.base_prefix)" OUTPUT_VARIABLE IN_VENV)
  string(STRIP ${IN_VENV} IN_VENV)
  # Maybe conda or mamba or ...
  execute_process(COMMAND ${Python_EXECUTABLE} -c "import sys, os; print(os.path.exists(os.path.join(sys.prefix, 'conda-meta')))" OUTPUT_VARIABLE IN_CONDA)
  string(STRIP ${IN_CONDA} IN_CONDA)

  if(SICONOS_SYSTEM_WIDE_INSTALL)
    # Keep default CMAKE_INSTALL_PREFIX (/usr/local ...)
    # Get python install path

    if(IN_VENV OR IN_CONDA)
      message(WARNING "You asked for a system wide install while cmake was run in a Python venv or conda environment. \n
Since make/pip install must be run as root, possibly outside this env., this might lead to very unexpected results.")
    endif()
    execute_process(COMMAND ${Python_EXECUTABLE} -c
      "import site; print(site.getsitepackages()[0])" OUTPUT_VARIABLE PY_INSTALL_DIR)
  else()
    if(IN_VENV)
      execute_process(COMMAND ${Python_EXECUTABLE} -c "import sys;print(sys.path[-1])" OUTPUT_VARIABLE PY_INSTALL_DIR)
    else()
      if(IN_CONDA AND DEFINED ENV{CONDA_PREFIX})
	execute_process(COMMAND ${Python_EXECUTABLE} -c "import sys;print(sys.path[-1])" OUTPUT_VARIABLE PY_INSTALL_DIR)
      else()
	execute_process(COMMAND ${Python_EXECUTABLE} -m site --user-base OUTPUT_VARIABLE PY_INSTALL_DIR)
      endif()
    endif()

    if(SICONOS_CUSTOM_INSTALL)
      set(CMAKE_INSTALL_PREFIX ${SICONOS_CUSTOM_INSTALL} CACHE PATH "Install root directory." FORCE)
      if(ISOLATED_INSTALL)
	list(APPEND PIP_INSTALL_OPTIONS_LOCAL --prefix=${SICONOS_CUSTOM_INSTALL})
	# Change python install path to custom
	set(PY_INSTALL_DIR ${SICONOS_CUSTOM_INSTALL}/lib/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/site-packages)
      else()
	# keep PY_INSTALL_DIR computed during check python env
	# There, we need to append --user to pip command
	list(APPEND PIP_INSTALL_OPTIONS_LOCAL --user)
      endif()
    else() # Default for any other cases
      # Are we in some kind of virtual env?
      if(IN_VENV)
	if(DEFINED ENV{VIRTUAL_ENV})
	  set(CMAKE_INSTALL_PREFIX $ENV{VIRTUAL_ENV} CACHE PATH "User install default directory" FORCE)
	else()
	  set(CMAKE_INSTALL_PREFIX $ENV{HOME}/.siconos CACHE PATH "User install default directory" FORCE)
	endif()
      else()
	if(IN_CONDA AND DEFINED ENV{CONDA_PREFIX})
	  set(CMAKE_INSTALL_PREFIX $ENV{CONDA_PREFIX} CACHE PATH "User install default directory" FORCE)
	else()
	  # There, we need to append --user to pip command
	  list(APPEND PIP_INSTALL_OPTIONS_LOCAL --user)
	  set(CMAKE_INSTALL_PREFIX $ENV{HOME}/.siconos CACHE PATH "User install default directory" FORCE)
	endif()
      endif()
    endif()
  endif()

  # Move vars SICONOS_PYTHON_INSTALL_DIR and PIP_INSTALL_OPTIONS to cache
  string(STRIP "${PY_INSTALL_DIR}" PY_INSTALL_DIR)
  set(SICONOS_PYTHON_INSTALL_DIR ${PY_INSTALL_DIR} CACHE PATH "Install directory for python bindings." FORCE)
  # Warning: this SICONOS_PYTHON_INSTALL_DIR will be used during swig setup to choose the place where dynamic libraries generated
  # by swig and required by python will be installed.
  list(APPEND PIP_INSTALL_OPTIONS_LOCAL --no-build-isolation) # To allow siconos install when network is not available
  set(PIP_INSTALL_OPTIONS ${PIP_INSTALL_OPTIONS_LOCAL} CACHE STRING "Options passed to pip during installation." FORCE)
endfunction()
