#======================================
#
# cmake setup for numerics for windows.
#======================================

# for LAPACKE, see http://icl.cs.utk.edu/lapack-for-windows/lapack/
IF(MSVC)
  APPEND_C_FLAGS("/DHAVE_LAPACK_CONFIG_H")
  APPEND_C_FLAGS("/DLAPACK_COMPLEX_STRUCTURE")
  APPEND_C_FLAGS("/DADD_")
  APPEND_CXX_FLAGS("/DHAVE_LAPACK_CONFIG_H")
  APPEND_CXX_FLAGS("/DLAPACK_COMPLEX_STRUCTURE")
  APPEND_CXX_FLAGS("/DADD_")
ENDIF(MSVC)
