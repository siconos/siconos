if(NOT WITH_PYTHON_WRAPPER)
  return()
endif()

# ------------- Tests -------------
if(WITH_TESTING)
  # About tox: https://tox.wiki/en/latest/
  # See tox.ini and pyproject.toml files.
  # What for? We can run the test suite with "tox ./wrap in the build dir."
  configure_file(tox.ini.in tox.ini)
  configure_file(conftest.py.in conftest.py) # to set data dir or other useful pytest fixtures.
  configure_file(requirements4tests.txt requirements4tests.txt)

endif()

# ====== Generates python files required for build/install process ======

# -- setuptools config --
# doc: https://packaging.python.org/en/latest/tutorials/packaging-projects/
# Required files: setup.cfg, MANIFEST.in, pyproject.toml.

configure_file(setup.cfg.in setup.cfg)
# A trick to add so files to the python package dependencies, so
# that they are installed by pip.
# This is a quite a mess to handle python path. See https://www.python.org/dev/peps/pep-0668/
file(GLOB libfiles CONFIGURE_DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/libs/*_lib)
foreach(fname IN LISTS libfiles)
  file(READ ${fname} VAR)
  file(APPEND  ${CMAKE_CURRENT_BINARY_DIR}/setup.cfg "\t${VAR}\n")
endforeach()

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

