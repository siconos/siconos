# --- Python ---

include(FindPythonInterp)
include(FindPythonLibs)
include(FindPythonModule)
if(NOT PYTHONLIBS_VERSION_STRING VERSION_EQUAL PYTHON_VERSION_STRING)
  display(PYTHONLIBS_VERSION_STRING)
  display(PYTHON_VERSION_STRING)
  message(FATAL_ERROR "Python library and executable versions do not match.")
endif()

find_python_module(numpy REQUIRED)
#Install dir for python (default = --user)
execute_process(
  COMMAND ${PYTHON_EXECUTABLE} -c "import site ; print site.USER_SITE"
  OUTPUT_VARIABLE ${PROJECT_NAME}_INSTALL_DIR)
string(STRIP ${${PROJECT_NAME}_INSTALL_DIR} ${PROJECT_NAME}_INSTALL_DIR)
set(CMAKE_INSTALL_PREFIX ${${PROJECT_NAME}_INSTALL_DIR}/${PROJECT_NAME})

# ============= Generates setup.py =============
# The file setup.py will be generated from setup.py.in.
if(EXISTS ${CMAKE_SOURCE_DIR}/setup.py.in)
  message(STATUS "Generate setup.py file ...")
  configure_file(setup.py.in setup.py)
endif()

set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES ${CMAKE_BINARY_DIR}/build)

add_custom_target(python-build ALL COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/setup.py build config_fc --f90exec=${CMAKE_Fortran_COMPILER}
   WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMENT "build parmepy package")

# To install python package AND parmes library and modules
if(WITH_LIB_FORTRAN)
  add_custom_target(python-install COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/setup.py install --user config_fc --f90exec=${CMAKE_Fortran_COMPILER}
    COMMAND make install
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMENT "build/install parmepy package")
else()
  add_custom_target(python-install COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/setup.py install --user config_fc --f90exec=${CMAKE_Fortran_COMPILER}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMENT "build/install parmepy package")
endif()

# Target to remove parmepy from install-path.
add_custom_target(python-cleaninstall COMMAND rm -rf ${${PROJECT_NAME}_INSTALL_DIR}/${PROJECT_NAME}*
  COMMENT "remove parmepy package and its dependencies")

# Target to clean sources (remove .pyc files) and build dir.
file(GLOB_RECURSE PYCFILES "${CMAKE_SOURCE_DIR}/*.pyc"})
add_custom_target(pyclean COMMAND rm -f ${PYCFILES}
  COMMAND make clean
  COMMAND rm -rf ${CMAKE_BINARY_DIR}/build  ${CMAKE_BINARY_DIR}/ParmesDoc
  COMMAND rm ${CMAKE_SOURCE_DIR}/parmepy/__init__.py
  COMMENT "clean parmepy sources and build.")

add_dependencies(python-build ${PARMES_LIBRARY_NAME})
add_dependencies(python-install ${PARMES_LIBRARY_NAME})

