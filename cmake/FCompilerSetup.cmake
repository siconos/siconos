# ===========================================
#  Set some predefined compilation flags for
#  the Fortran compiler.
#
#  Usage :
#  include(FCompilerSetup)
# 
# ===========================================

include(FortranCInterface)
# Set module files directory (i.e. where .mod will be created)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Modules)

append_Fortran_FLAGS("-fPIC")
if(DEV_MODE)
  append_Fortran_FLAGS("-w") # gnu specific ...
endif()