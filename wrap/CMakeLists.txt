# ============================================================
# Prepare swig config to generate python bindings for siconos
# ============================================================
if(WITH_PYTHON_WRAPPER)
  # Name of the generated Python package
  set(SICONOS_PYTHON_PACKAGE ${PROJECT_NAME} CACHE INTERNAL "Name of the siconos python module.")
  
  # --------------- Python install setup ---------------
  # Directory of all config files for python bindings.
  # Set path for siconos-python installation
  include(PythonInstallSetup)
  set_python_install_path()
  
  #  ----------- swig -----------

  # -- Swig flags --
  if(SYSTEM_ARCHITECTURE)
    list(APPEND CMAKE_SWIG_FLAGS "-D__${SYSTEM_ARCHITECTURE}__")
  endif()

  if(PYTHON_VERSION_MAJOR VERSION_GREATER 2)
    list(APPEND CMAKE_SWIG_FLAGS "-py3")
  endif()

  list(REMOVE_DUPLICATES CMAKE_SWIG_FLAGS)
  set(CMAKE_SWIG_FLAGS "${CMAKE_SWIG_FLAGS}" CACHE INTERNAL "Swig flags")

  # -dirvtable      - Generate a pseudo virtual table for directors for faster dispatch
  # - without "dirprot", swig wrap public methods and only the 
  # protected methods needed to the interface to compile.
  # - with "dirprot" swig will attemp to wrap all the public and protected methods at once.
  set(SWIG_DEFS "-dirprot;-dirvtable;-Wextra")
  
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
  # i.e. .i files, common to all modules and submodules.
  set(SICONOS_SWIG_SRC_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/swig
    CACHE INTERNAL "Path to swig files common to all packages.")

  # -- Output dir --
  set(SICONOS_SWIG_ROOT_DIR ${CMAKE_CURRENT_BINARY_DIR}/${SICONOS_PYTHON_PACKAGE}
    CACHE INTERNAL "Root path for swig outputs.")

  # -- include dir --
  set(SICONOS_SWIG_INCLUDE_DIRS ${SICONOS_SWIG_SRC_DIRS} ${SICONOS_SWIG_ROOT_DIR} ${PYTHON_INCLUDE_DIRS}
    CACHE INTERNAL "Directories required for swig includes.")

  foreach(_dir ${SICONOS_SWIG_INCLUDE_DIRS})
    set(SWIG_DEFS "-I${_dir};${SWIG_DEFS}")
  endforeach()

  set(SWIG_DEFS "${SWIG_DEFS}" CACHE INTERNAL "Swig extra definitions.")

  include_directories(${SICONOS_SWIG_INCLUDE_DIRS})

  # ------------- Tests -------------
  if(WITH_TESTING)
    IF(CROSSCOMPILING_LINUX_TO_WINDOWS)
      SET(EMULATOR "wine")
      SET(DRIVE_LETTER "Z:")
    ELSE()
      SET(EMULATOR)
      SET(DRIVE_LETTER)
    ENDIF()
    
    find_python_module(pytest REQUIRED)
    # -- Runner for tests --
    set(TESTS_RUNNER ${CMAKE_CURRENT_BINARY_DIR}/runtests.py CACHE INTERNAL "driver for python tests")
    if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/runtests.py.in)
      message(STATUS "Generate runtests.py file ...")
      configure_file(runtests.py.in ${TESTS_RUNNER})
    endif()
    if(WITH_AGGRESSIVE_PYTHON_TESTS)
      set(pytest_opt "-s -v -pep8" CACHE INTERNAL "extra options for py.test")
    else()
      set(pytest_opt "-v" CACHE INTERNAL "extra options for py.test")
    endif()
    file(MAKE_DIRECTORY ${SICONOS_SWIG_ROOT_DIR}/tests)
    configure_file(tests_setup.py.in ${SICONOS_SWIG_ROOT_DIR}/tests_setup.py)
    configure_file(set_tests_env.sh.in set_tests_env.sh)
    configure_file(utils/__init__.py ${SICONOS_SWIG_ROOT_DIR}/tests/__init__.py)
    add_custom_target(set_envtest numerics
      COMMAND sh ${CMAKE_CURRENT_BINARY_DIR}/set_tests_env.sh
      COMMENT "set environment variables for python tests ...")
    # todo if we need to add some 'global' tests for siconos (not component dependent)
    # create a wrap/swig/tests dir and ...
    #build_python_tests()
  endif()

  # ====== Generates python files required for build/install process ======
  
  set(SICONOS_PYTHON_MODULES)
  foreach(comp ${COMPONENTS})
    if(WITH_${comp}_PYTHON_WRAPPER)
      set(SICONOS_PYTHON_MODULES "'${comp}', ${SICONOS_PYTHON_MODULES}")
    endif()
  endforeach()
  # setup.py --> CMAKE_BINARY_DIR
  if(EXISTS  ${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in)
    message(STATUS "Generate setup.py file ...")
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in setup.py)
  endif()

  # The file __init__.py will be generated from __init__.py.in.
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py.in)
    message(STATUS "Generate __init__.py file ...")
    if(EXISTS ${SICONOS_SWIG_ROOT_DIR}/__init__.py)
      file(REMOVE ${SICONOS_SWIG_ROOT_DIR}/__init__.py)
    endif()
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/__init__.py.in ${SICONOS_SWIG_ROOT_DIR}/__init__.py)
  endif()

  # Some python files, that do not depend on siconos
  # but may be useful ...
  file(GLOB pyfiles ${CMAKE_CURRENT_SOURCE_DIR}/utils/[^_]*.py)
  foreach(_file ${pyfiles})
    get_filename_component(_filename ${_file} NAME)
    configure_file(utils/${_filename} ${SICONOS_SWIG_ROOT_DIR}/${_filename} COPYONLY)
  endforeach()

  # ====== Create (and setup) build target ======

  add_custom_target(python-build
    COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/setup.py build
    COMMENT "Build ${SICONOS_PYTHON_PACKAGE} package")

  IF(UNIX AND NOT APPLE)
    LIST(APPEND python_install_options "--root=$(DESTDIR)/")
  ENDIF(UNIX AND NOT APPLE)
  add_custom_target(python-install
    COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/setup.py install ${python_install_options}
    VERBATIM
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMENT "build/install siconos package")
  install(CODE "execute_process(COMMAND ${CMAKE_BUILD_TOOL} python-install WORKING_DIRECTORY \"${CMAKE_CURRENT_BINARY_DIR}\")")

  # ====== Create a Target to clean sources (remove .pyc files) and build dir ======

  file(GLOB_RECURSE PYCFILES "${CMAKE_SOURCE_DIR}/*.pyc")
  add_custom_target(pyclean COMMAND rm -f ${PYCFILES}
    COMMAND make clean
    COMMAND rm -rf ${CMAKE_BINARY_DIR}/build
    COMMENT "Clean sources and build.")

  # ====== conf/copy file for docstrings generation from doxygen ====
  if(WITH_DOXY2SWIG)
    configure_file(doxy2swig.py ${CMAKE_BINARY_DIR}/share/doxy2swig.py COPYONLY)
  endif()

endif()
