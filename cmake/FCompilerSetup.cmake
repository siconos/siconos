# ===========================================
#  Set some predefined compilation flags for
#  the Fortran compiler.
#
#  Usage :
#
#   project(... Fortran ...)
#   include(FCompilerSetup) 
# 
# ===========================================

# Fortran needs C ...
enable_language(C)
include(FortranCInterface)

if(WITH_CXX)
  FortranCInterface_VERIFY(CXX QUIET)
else()
  FortranCInterface_VERIFY(QUIET)
endif()


if(NOT FortranCInterface_VERIFIED_C)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/DefaultMangle.h.in ${CMAKE_CURRENT_BINARY_DIR}/FCMangle.h)
else()
  FortranCInterface_HEADER(${CMAKE_CURRENT_BINARY_DIR}/FCMangle.h
    MACRO_NAMESPACE "myF2C" 
    SYMBOL_NAMESPACE "myF2C"
    )
endif()

# Set module files directory (i.e. where .mod will be created)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Modules)
