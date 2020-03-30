#======================================
#
# cmake setup for numerics for windows.
#======================================

# for LAPACKE, see http://icl.cs.utk.edu/lapack-for-windows/lapack/
if(MSVC)
  target_compile_definitions(numerics PRIVATE /DHAVE_LAPACK_CONFIG_H)
  target_compile_definitions(numerics PRIVATE /DLAPACK_COMPLEX_STRUCTURE)
  target_compile_definitions(numerics PRIVATE /DADD_)
endif()
