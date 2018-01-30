# ================================================================
# All the default values for siconos cmake parameters
#
# Usage:
# cmake path-to-sources
#  --> to keep default value
# 
# cmake path-to-sources -DWITH_PYTHON_WRAPPER=ON
#  --> to enable (ON), or disable (OFF) the concerned option.
#
# For details about all these options check siconos install guide.
# ================================================================

# --------- User-defined options ---------
# Use cmake -DOPTION_NAME=some-value ... to modify default value.
option(WITH_DOCUMENTATION "Build Documentation. Default = OFF" OFF)
option(WITH_PYTHON_WRAPPER "Build python bindings using swig. Default = ON" ON)
option(WITH_DOXYGEN_WARNINGS "Explore doxygen warnings. Default = OFF" OFF)
option(WITH_DOXY2SWIG "Build swig docstrings from doxygen xml output. Default = OFF." OFF)
option(WITH_SYSTEM_INFO "Verbose mode to get some system/arch details. Default = OFF." OFF)
option(WITH_TESTING "Enable 'make test' target" OFF)
option(WITH_GIT "Consider sources are under GIT" OFF)
option(WITH_SERIALIZATION "Compilation of serialization functions. Default = OFF" OFF)
option(WITH_GENERATION "Generation of serialization functions with doxygen XML. Default = OFF" OFF)
option(WITH_CXX "Enable CXX compiler for numerics. Default = ON" ON)
option(WITH_FORTRAN "Enable Fortran compiler. Default = ON" ON)
option(WITH_UNSTABLE "Enable this to include all 'unstable' sources. Default=OFF" OFF)
option(WITH_UNSTABLE_TEST "Enable this to include all 'unstable' test. Default=OFF" OFF)
option(BUILD_SHARED_LIBS "Building of shared libraries. Default = ON" ON)
option(DEV_MODE "Compilation flags setup for developers. Default = OFF" OFF)
option(DEV_MODE_STRICT "Compilation flags setup for developers (extra strict, conversion warnings). Default = OFF" OFF)
option(WITH_BULLET "compilation with Bullet Bindings. Default = OFF" OFF)
option(WITH_OCC "compilation with OpenCascade Bindings. Default = OFF" OFF)
option(WITH_MUMPS "Compilation with the MUMPS solver. Default = OFF" OFF)
option(WITH_UMFPACK "Compilation with the UMFPACK solver. Default = OFF" OFF)
option(WITH_SUPERLU "Compilation with the SuperLU solver. Default = OFF" OFF)
option(WITH_SUPERLU_MT "Compilation with the SuperLU solver, multithreaded version. Default = OFF" OFF)
option(WITH_FCLIB "link with fclib when this mode is enable. Default = OFF" OFF)
option(WITH_FREECAD "Use FreeCAD. Default = OFF" OFF)
option(WITH_MECHANISMS "Generation of bindings for Mechanisms toolbox (required OCE). Default = OFF" OFF)
option(WITH_RENDERER "Install OCC renderer. Default = OFF" OFF)
option(WITH_SYSTEM_SUITESPARSE "Use SuiteSparse installed on the system instead of built-in CXSparse library. Default = ON" ON)
option(WITH_XML "Enable xml files i/o. Default = ON" ON)
option(WITH_DOCKER "Build inside a docker container. Default = OFF" OFF)
option(FORCE_SKIP_RPATH "Do not build shared libraries with rpath. Useful only for packaging. Default = OFF" OFF)
option(NO_RUNTIME_BUILD_DEP "Do not check for runtime dependencies. Useful only for packaging. Default = OFF" OFF)

# Set python install mode:
# - user --> behave as 'python setup.py install --user'
# - standard --> install in python site-package (ie behave as python setup.py install)
# - prefix --> install in python CMAKE_INSTALL_PREFIX (ie behave as python setup.py install --prefix=CMAKE_INSTALL_PREFIX)
IF(UNIX)
  # on unix, there is no reason to use the standard option. By default, CMAKE_INSTALL_PREFIX is set to /usr/local and therefore,
  # the python packages should be installed in /usr/local/...
  set(siconos_python_install "prefix" CACHE STRING "Install mode for siconos python package")
ELSE(UNIX)
  set(siconos_python_install "standard" CACHE STRING "Install mode for siconos python package")
ENDIF(UNIX)

# If OFF, headers from libraries in externals will not be installed.
option(INSTALL_EXTERNAL_HEADERS
  "Whether or not headers for external libraries should be installed. Default=OFF" OFF)

# If ON, internal headers will not be installed.
option(INSTALL_INTERNAL_HEADERS
  "Whether or not headers for internal definitions should be installed. Default=OFF" OFF)

# List of components to build and installed
# List of siconos component to be installed
# complete list = externals numerics kernel control mechanics io
set(COMPONENTS externals numerics kernel control mechanics io CACHE INTERNAL "List of siconos components to build and install")
