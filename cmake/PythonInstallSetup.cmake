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
#[=======================================================================[.rst:
PythonInstallSetup
------------------
 --- Function to compute the install path/options ---

 Usage :
 set_install_options()

 - Sets PIP_INSTALL_OPTIONS (cache)
 Reads siconos_python_install option (from default_user_options)
 and set python install dir according to its value
 
 - Sets SICONOS_PYTHON_INSTALL_DIR (cache) with the full path where
 Siconos python packages will be installed.

 
 Summary :
 cmake path-to-your-sources -Dpython_install_dir=standard
 make install

 * install_dir = 'user'
   ---> install python packages with pip install --user 

 * install_dir = 'prefix'
   ---> install python packages with pip install --prefix=$CMAKE_INSTALL_PREFIX

 * else
   ---> install python packages with pip install ...
#]=======================================================================]


function(set_python_install_path)
  # set(PIP_INSTALL_OPTIONS "--record;${CMAKE_BINARY_DIR}/python_install_manifest.txt")
  set(PIP_INSTALL_OPTIONS)
  execute_process(COMMAND ${PYTHON_EXECUTABLE} -c
    "import sys; print('%d.%d'%(sys.version_info.major,sys.version_info.minor))"
    OUTPUT_VARIABLE PY_VERSION)
  string(STRIP ${PY_VERSION} PY_VERSION)

  if(siconos_python_install STREQUAL "user")
    # --- Case 1 : siconos_python_install=user ---
    # Will add "--user" option to pip install
    # In that case, we need to find the user path. It depends on the operation system
    # and on which python is used (virtualenv ...)
    # This path is only required by cmake to install generated libraries (swig)
    # (see swig_python_tools.cmake)
    # Find install path for --user (site.USER_SITE)

    # -- PY_INSTALL_DIR 
    execute_process(COMMAND ${PYTHON_EXECUTABLE} -c
      "import site; print(site.ENABLE_USER_SITE)" OUTPUT_VARIABLE ENABLE_USER_SITE)
    string(STRIP ${ENABLE_USER} ENABLE_USER)

    if(ENABLE_USER_SITE)
      execute_process(COMMAND ${PYTHON_EXECUTABLE} -c
	"import site; print(site.USER_SITE)" OUTPUT_VARIABLE PY_INSTALL_DIR)

      # -- pip options
      list(APPEND PIP_INSTALL_OPTIONS --user)
    else()
      message(FATAL_ERROR "pip (python install) can't use --user option. Please change your install configuration setup (siconos_python_install var).")
    endif()

  elseif(siconos_python_install STREQUAL prefix) # User explicitely asked to use CMAKE_INSTALL_PREFIX
    # -- pip options
    list(APPEND PIP_INSTALL_OPTIONS --prefix=${CMAKE_INSTALL_PREFIX})

    # -- PY_INSTALL_DIR 
    execute_process(COMMAND ${PYTHON_EXECUTABLE} -c
      "import site; print('/lib'+site.getsitepackages()[0].split('lib')[-1])" OUTPUT_VARIABLE PY_INSTALL_DIR)

    set(PY_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}${PY_INSTALL_DIR})

  else()  # Default behavior: let pip choose where things will be installed.
    # Anyway, we need PY_INSTALL_DIR for swig/cmake install process.
    execute_process(COMMAND ${PYTHON_EXECUTABLE} -c
      "import site; print(site.getsitepackages()[0])" OUTPUT_VARIABLE PY_INSTALL_DIR)
    
  endif()
  # Move vars to cache.
  string(STRIP ${PY_INSTALL_DIR} PY_INSTALL_DIR)
  set(SICONOS_PYTHON_INSTALL_DIR ${PY_INSTALL_DIR} CACHE PATH "Install directory for python bindings.")
  #string(STRIP ${python_install_options} python_install_options)
  set(PIP_INSTALL_OPTIONS ${PIP_INSTALL_OPTIONS} CACHE STRING "Options passed to pip during installation.")
endfunction()

