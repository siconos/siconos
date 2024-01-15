# =========================================================
# Advanced options used by cmake to configure Siconos
#
# The default values are sufficient for most applications.
#
# Thus, the following options are mostly for
# developers or advanced users and so hidden here.
#
# They could be easyly updated by using
#
# cmake -DOPTION_NAME=NEW_VALUE ...
#
# =========================================================
include(CMakeDependentOption)

option(WITH_SERIALIZATION "Compilation of serialization functions" OFF)
option(WITH_GENERATION "Generation of serialization functions with doxygen XML" OFF)

# --- Build/compiling options ---
option(WITH_CXX "Enable CXX compiler for numerics" ON)
option(WITH_FORTRAN "Enable Fortran compiler" ON)

option(FORCE_SKIP_RPATH "Do not build shared libraries with rpath. Useful only for packaging" OFF)
option(NO_RUNTIME_BUILD_DEP "Do not check for runtime dependencies. Useful only for packaging" OFF)
option(WITH_UNSTABLE_TEST "Enable this to include all 'unstable' test. Default=OFF" OFF)
option(BUILD_SHARED_LIBS "Building of shared libraries" ON)
option(WITH_SYSTEM_INFO "Verbose mode to get some system/arch details." OFF)

# --- Documentation setup ---
cmake_dependent_option(WITH_DOXY2SWIG  "Build swig docstrings from doxygen xml output if doc is ON." ON "WITH_DOCUMENTATION" OFF)
cmake_dependent_option(WITH_DOXYGEN_WARNINGS "Explore doxygen warnings if doc is ON" ON "WITH_DOCUMENTATION" OFF)

# --- List of external libraries/dependencies to be searched (or not) ---
option(WITH_MA57 "Compilation with the MA57 solver (License HSL)" OFF)
option(WITH_FREECAD "Use FreeCAD" OFF)
option(WITH_RENDERER "Install OCC renderer" OFF)
option(WITH_SYSTEM_SUITESPARSE "Use SuiteSparse installed on the system instead of built-in CXSparse library" ON)
option(WITH_XML "Enable xml files i/o" OFF)

# If OFF, headers from libraries in externals will not be installed.
option(INSTALL_EXTERNAL_HEADERS
  "Whether or not headers for external libraries should be installed" OFF)

# If ON, internal headers will not be installed.
option(INSTALL_INTERNAL_HEADERS
  "Whether or not headers for internal definitions should be installed" OFF)
