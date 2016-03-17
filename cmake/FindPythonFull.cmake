#.rst:
# FindPythonFull
# --------------
#
# Find python interpreter and all required libraries and headers.
#
# The default cmake find_package(PythonLibs) process does not work for us.
#
# Usage:
# find_package(PythonFull)
#
# This call will set the following variables :
# ::
#
#   PYTHON_FOUND           - True if Python executable, libraries and header were found
#   PYTHON_EXECUTABLE          - python interpreter (full path)
#   PYTHON_VERSION_STRING      - Python version found e.g. 2.5.2
#   PYTHON_LIBRARIES           - full path to the python library
#   PYTHON_INCLUDE_DIRS        - full path to Python.h
#   PYTHONLIBS_VERSION_STRING  - version of the Python libs found
#
# By default, we search for the current active python version first.
# If you need another version, use -DPYTHON_EXECUTABLE=full-path-to-python-exe
# during cmake call.
#

set(PYTHON_FOUND FALSE)

# Does nothing if vars are already in cache
if(EXISTS "${PYTHON_INCLUDE_DIRS}" AND EXISTS "${PYTHON_LIBRARY}" AND EXISTS "${PYTHON_SITE_PACKAGES_DIR}")
  set(PYTHON_FOUND TRUE)
else()
  set(PYTHON_FOUND FALSE)
  # --- Find python interpreter
  find_package(PythonInterp)

  # --- Use distutils to explore python configuration corresponding to
  # the python executable found.
  find_file(_findpython explore_python_config.py PATHS ${CMAKE_MODULE_PATH})

  execute_process(
    COMMAND ${PYTHON_EXECUTABLE} ${_findpython}
    OUTPUT_VARIABLE python_config
    )

  # --- Post-process distutils results
  if(python_config)
    string(REGEX REPLACE ".*exec_prefix:([^\n]+).*$" "\\1" PYTHON_PREFIX ${python_config})
    string(REGEX REPLACE ".*\nversion:([^\n]+).*$" "\\1" PYTHON_VERSION ${python_config})
    string(REGEX REPLACE ".*\npy_inc_dir:([^\n]+).*$" "\\1" PYTHON_INCLUDE_DIRS ${python_config})
    string(REGEX REPLACE ".*\nsite_packages_dir:([^\n]+).*$" "\\1" PYTHON_SITE_PACKAGES_DIR ${python_config})
    string(REGEX REPLACE "([0-9]+).([0-9]+)" "\\1\\2" PYTHON_VERSION_NO_DOTS ${PYTHON_VERSION})
    if(WIN32)
      string(REPLACE "\\" "/" PYTHON_SITE_PACKAGES_DIR ${PYTHON_SITE_PACKAGES_DIR})
      string(REPLACE "\\" "/" PYTHON_PREFIX ${PYTHON_PREFIX})
    endif(WIN32)

    set(PY_VARIANTS "mu;m;u")
    set(PYLIB_NAMES python${PYTHON_VERSION_NO_DOTS})
    set(PYLIB_HINTS ${PYTHON_PREFIX})
    foreach(_PY_V ${PY_VARIANTS})
      set(_PY_N  python${PYTHON_VERSION}${_PY_V})
      list(APPEND PYLIB_NAMES ${_PY_N})
      list(APPEND PYLIB_HINTS ${PYTHON_PREFIX}/lib/${_PY_N}/config) # I'm not sure this one is used anywhere, but it doesn't hurt to add it
      list(APPEND PYLIB_HINTS ${PYTHON_PREFIX}/lib/python${PYTHON_VERSION}/config-${PYTHON_VERSION}${_PY_V}-${CMAKE_LIBRARY_ARCHITECTURE})
    endforeach(_PY_V ${PY_VARIANTS})
    list(APPEND PYLIB_NAMES python${PYTHON_VERSION})
    list(APPEND PYLIB_HINTS ${PYTHON_PREFIX}/lib/python${PYTHON_VERSION}/config)
    list(APPEND PYLIB_HINTS ${PYTHON_PREFIX}/lib/python${PYTHON_VERSION}/config-${CMAKE_LIBRARY_ARCHITECTURE})

    # --- Search python library corresponding to python exec.
    find_library(PYTHON_LIBRARY
      NAMES ${PYLIB_NAMES}
      NO_DEFAULT_PATH
      HINTS ${PYLIB_HINTS}
      PATH_SUFFIXES lib libs
      )

    set(PYTHON_LIBRARIES ${PYTHON_LIBRARY} CACHE FILEPATH "Python libraries" FORCE)

    set(PYTHON_INCLUDE_DIRS ${PYTHON_INCLUDE_DIRS} CACHE FILEPATH "Path to Python.h" FORCE)

    # --- Extract python library version for further checks.
    if(PYTHON_INCLUDE_DIRS AND EXISTS "${PYTHON_INCLUDE_DIRS}/patchlevel.h")
      file(STRINGS "${PYTHON_INCLUDE_DIRS}/patchlevel.h" python_version_str
        REGEX "^#define[ \t]+PY_VERSION[ \t]+\"[^\"]+\"")
      string(REGEX REPLACE "^#define[ \t]+PY_VERSION[ \t]+\"([^\"]+)\".*" "\\1"
        PYTHONLIBS_VERSION_STRING "${python_version_str}")
      unset(python_version_str)
    endif()
    
  endif()

  unset(PYTHON_FOUND)
  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(Python
    REQUIRED_VARS PYTHON_LIBRARIES PYTHON_INCLUDE_DIRS PYTHON_EXECUTABLE
    VERSION_VAR PYTHONLIBS_VERSION_STRING)
  if(PYTHON_FOUND)
    set(PYTHONFULL_FOUND TRUE)
    if(NOT PythonFull_FIND_QUIETLY)
      message("-- Found Python executable: ${PYTHON_EXECUTABLE}")
      message("-- Found Python library: ${PYTHON_LIBRARIES}")
      message("-- Python version is : ${PYTHON_VERSION_STRING}")
      message("-- Python include dir is : ${PYTHON_INCLUDE_DIRS}")
      message("-- Python Site package dir is : ${PYTHON_SITE_PACKAGES_DIR}\n")
    endif()
  else()
    if(PythonFull_FIND_REQUIRED)
      message(FATAL_ERROR "Could not find Python")
    endif()
  endif()
  
endif()

if(NOT PYTHONLIBS_VERSION_STRING VERSION_EQUAL PYTHON_VERSION_STRING)
  print_var(PYTHONLIBS_VERSION_STRING)
  print_var(PYTHON_VERSION_STRING)
  message(FATAL_ERROR "Python library and executable versions do not match. Please check your python installation.")
endif()
