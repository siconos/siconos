if(NOT WITH_PYB11_WRAPPER)
  return()
endif()

message("------------------------------------------\n")
message("-- Start set up for Python Siconos API ...\n")

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

configure_file(pyproject.toml pyproject.toml)

if(EXISTS  ${CMAKE_CURRENT_SOURCE_DIR}/MANIFEST.in)
  configure_file(MANIFEST.in MANIFEST.in)
endif()

# packaging
configure_file(${CMAKE_SOURCE_DIR}/COPYING LICENSE)
configure_file(${CMAKE_SOURCE_DIR}/README.md README.md)

# -- siconos package --
message(STATUS "Generate __init__.py file ...")

configure_file(__init__.py.in ${CMAKE_CURRENT_BINARY_DIR}/siconos/__init__.py)

# ====== Create a Target to clean sources (remove .pyc files) and build dir ======
file(GLOB_RECURSE PYCFILES "${CMAKE_SOURCE_DIR}/*.pyc")
add_custom_target(pyclean COMMAND rm -f ${PYCFILES}
  COMMAND make clean
  COMMAND rm -rf ${CMAKE_BINARY_DIR}/build
  COMMENT "Clean sources and build.")

message("-- End set up for Python Siconos API\n")
message("------------------------------------------\n")
