
option(INSTALL_PYTHON_SYMLINKS "Install Python .py files as symlinks" OFF)

# --- Function to compute the install path/options ---
#
# Usage :
# set_install_options()
#
# Reads siconos_python_install option (from default_user_options)
# and set python install dir according to its value.
#
# Summary :
# cmake path-to-your-sources -Dpython_install_dir=standard
# make install
#
# ---> install in 'system' python site-package
#
# cmake path-to-your-sources -Dpython_install_dir=user
# make install
#
# ---> install in USER_SITE (no virtualenv case)
# ---> install in site-packages of your virtualenv
#
# cmake path-to-your-sources -Dpython_install_dir=prefix -DCMAKE_INSTALL_PREFIX=/some/install/path
# make install
#
# ---> install in CMAKE_INSTALL_PREFIX
#
# If /some/install/path is not a standard path of your system,
# you'll probably need something like :
# export PYTHONPATH=${PYTHONPATH}:/some/install/path
#
#
function(set_python_install_path)
  set(python_install_options "--record;${CMAKE_BINARY_DIR}/python_install_manifest.txt")
  if(siconos_python_install STREQUAL "user")
    # --- Case 1 : siconos_python_install=user ---
    # In that case, we need to find the user path. It depends on the operation system
    # and on which python is used (virtualenv or not)
    # First, we need to check if '--user' option works in the current environment.
    execute_process(COMMAND ${PYTHON_EXECUTABLE} -c
      "import site; print(site.ENABLE_USER_SITE)" OUTPUT_VARIABLE ENABLE_USER)
    string(STRIP ${ENABLE_USER} ENABLE_USER)
    
    if(ENABLE_USER) # --user works ...
      # Find install path for --user (site.USER_SITE)
      execute_process(COMMAND ${PYTHON_EXECUTABLE} -c
       "import site; print(site.USER_BASE)" OUTPUT_VARIABLE USER_BASE)
      list(APPEND python_install_options --prefix=${USER_BASE})
      # Get python user site and install path = USER_SITE + project_name
      set(PYTHON_COMMAND_GET_INSTALL_DIR
       "import site, os, sys ; print(os.path.join(site.USER_BASE, os.path.join(\"lib\", os.path.join(\"python\" + str(sys.version_info.major) + '.' + str(sys.version_info.minor),
 \"site-packages\"))))")
    execute_process(
      COMMAND ${PYTHON_EXECUTABLE} -c "${PYTHON_COMMAND_GET_INSTALL_DIR}"
      OUTPUT_VARIABLE PY_INSTALL_DIR)

    else()
      # user site not included in the path,
      # which probably means that python is run using virtualenv
      # Command to find 'global' site-packages
      # default path will probably be ok --> no options
      set(GET_SITE_PACKAGE
       "from distutils.sysconfig import get_python_lib; print(get_python_lib())")
      execute_process(COMMAND ${PYTHON_EXECUTABLE} -c
       "${GET_SITE_PACKAGE}" OUTPUT_VARIABLE PY_INSTALL_DIR)
    endif()
    # Set the SICONOS_PYTHON_INSTALL_DIR to the proper path

  elseif(siconos_python_install STREQUAL prefix)
    # Case 2 : siconos_python_install=prefix
    # we use CMAKE_INSTALL_PREFIX as the path for python install
    list(APPEND python_install_options --prefix=${CMAKE_INSTALL_PREFIX})
      set(GET_SITE_PACKAGE
       "from distutils.sysconfig import get_python_lib; print(get_python_lib(prefix='${CMAKE_INSTALL_PREFIX}'))")
      execute_process(COMMAND ${PYTHON_EXECUTABLE} -c
       "${GET_SITE_PACKAGE}" OUTPUT_VARIABLE PY_INSTALL_DIR)
  else()
    # Default case : siconos_python_install=standard
    #set(PYTHON_COMMAND_GET_INSTALL_DIR   
    #  "import site; print(site.getsitepackages()[0])")
    # --> this does not work properly: the order in resulting
    # list depends on the OS, the python version ...
    # Set the SICONOS_PYTHON_INSTALL_DIR to the proper path
    configure_file(fake_setup.py tmp/setup.py)
    configure_file(fake/__init__.py tmp/fake/__init__.py)
    configure_file(find_python_install.py tmp/find_python_install.py)
    execute_process(
      COMMAND ${PYTHON_EXECUTABLE} find_python_install.py
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/tmp/
      OUTPUT_VARIABLE PY_INSTALL_DIR)
  endif()
  string(STRIP ${PY_INSTALL_DIR} PY_INSTALL_DIR)
  set(SICONOS_PYTHON_INSTALL_DIR ${PY_INSTALL_DIR}
    CACHE PATH "Install directory for python bindings." FORCE)
  #string(STRIP ${python_install_options} python_install_options)
  set(python_install_options ${python_install_options} CACHE INTERNAL "")
endfunction()

