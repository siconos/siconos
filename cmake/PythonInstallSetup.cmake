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
 Reads sicopy_install_mode option (from default_user_options)
 and set python install dir according to its value
 
 - Sets SICONOS_PYTHON_INSTALL_DIR (cache) with the full path where
 Siconos python packages will be installed.

 
 Summary :
 cmake path-to-your-sources -Dsiconos_python_install=standard
 make install

 * siconos_python_install= 'user'
   ---> install python packages with pip install --user 

 * siconos_python_install = 'isolated'
   ---> install python packages with pip install --prefix=$CMAKE_INSTALL_PREFIX

 * else (standard)
   ---> install python packages with pip install ...
#]=======================================================================]
  
# Finding where pip or python will put things is a real nightmare, especially on Debian
# like systems where pip use site-packages, the distro dist-packages and so on.
# See https://www.python.org/dev/peps/pep-0668/
# Note : witth python 3.10, --prefix option of pip always add "local" to install path.
# We can not handle this in a simple way, so we just keep standard and user options.
function(set_python_install_path)
  # set(PIP_INSTALL_OPTIONS "--record;${CMAKE_BINARY_DIR}/python_install_manifest.txt")
  set(PIP_INSTALL_OPTIONS_LOCAL)

  if(sicopy_install_mode STREQUAL "user")
    # --- Case 1 : sicopy_install_mode=user ---
    # Will add "--user" option to pip install
    # In that case, we need to find the user path. It depends on the operation system
    # and on which python is used (virtualenv ...)
    # This path is only required by cmake to install generated libraries (swig)
    # (see swig_python_tools.cmake)
    # Find install path for --user (site.USER_SITE)

    if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
      # In that case, we can not allow "user" option for python install:
      # su/sudo run is usually required for default install path and 
      # it will break the installation process with some files in (e.g.) ~root/.local/lib/python/... and
      # some in user home/.local and so on.
      message(FATAL_ERROR "pip (python install) can't use --user option when using a default CMAKE_INSTALL_PREFIX which requires root privileges. Please either set CMAKE_INSTALL_PREFIX to a path on which you're authorized to write or use sicopy_install_mode 'standard'.")
    endif()
    
    # -- PY_INSTALL_DIR 
    execute_process(COMMAND ${Python_EXECUTABLE} -c
      "import site; print(site.ENABLE_USER_SITE)" OUTPUT_VARIABLE ENABLE_USER_SITE)
    string(STRIP ${ENABLE_USER_SITE} ENABLE_USER_SITE)

    if(ENABLE_USER_SITE)
      execute_process(COMMAND ${Python_EXECUTABLE} -c
	"import site; print(site.USER_SITE)" OUTPUT_VARIABLE PY_INSTALL_DIR)

      # -- pip options
      list(APPEND PIP_INSTALL_OPTIONS_LOCAL --user)
    else()
      message(FATAL_ERROR "pip (python install) can't use --user option. Please change your install configuration setup (sicopy_install_mode var).")
    endif()
  elseif(sicopy_install_mode STREQUAL "isolated")
    list(APPEND PIP_INSTALL_OPTIONS_LOCAL --prefix=${ISOLATED_INSTALL})
    set(PY_INSTALL_DIR ${ISOLATED_INSTALL}/lib/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/site-packages)
  elseif(sicopy_install_mode STREQUAL "standard")  # Default behavior: let pip choose where things will be installed.
    # Anyway, we need PY_INSTALL_DIR for swig/cmake install process.
    execute_process(COMMAND ${Python_EXECUTABLE} -c
      "import site; print(site.getsitepackages()[0])" OUTPUT_VARIABLE PY_INSTALL_DIR)

  else()
    message(FATAL_ERROR "Unknown python install mode ${sicopy_install_mode}. Choose between 'user', 'isolated' or 'standard'.")
  endif()
  # Move vars to cache.
  string(STRIP "${PY_INSTALL_DIR}" PY_INSTALL_DIR)
  set(SICONOS_PYTHON_INSTALL_DIR ${PY_INSTALL_DIR} CACHE PATH "Install directory for python bindings." FORCE)
  # Warning: this SICONOS_PYTHON_INSTALL_DIR will be used during swig setup to choose the place where dynamic libraries generated
  # by swig and required by python will be installed.
  list(APPEND PIP_INSTALL_OPTIONS_LOCAL --no-build-isolation) # To allow siconos install when network is not available
  set(PIP_INSTALL_OPTIONS ${PIP_INSTALL_OPTIONS_LOCAL} CACHE STRING "Options passed to pip during installation." FORCE)
endfunction()

