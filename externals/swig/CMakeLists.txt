#================================================================
# cmake utilities to build python-bindings for externals component
#================================================================

# wrap - Set extra dependencies
set(SWIG_MODULE_${COMPONENT}_EXTRA_DEPS
  ${SICONOS_SWIG_SRC_DIRS}/KernelTypes.i
  ${SICONOS_SWIG_SRC_DIRS}/solverOptions.i
  ${SICONOS_SWIG_SRC_DIRS}/start.i
  ${SICONOS_SWIG_SRC_DIRS}/end.i
  ${SICONOS_SWIG_SRC_DIRS}/numpy.i
  ${SICONOS_SWIG_SRC_DIRS}/pyRegister.i)

# List of modules to be built in the python bindings
set(${COMPONENT}_PYTHON_MODULES "./${COMPONENT}")

set(${COMPONENT}_SWIG_INCLUDE_DIRECTORIES)
set(${COMPONENT}_SWIG_DEFS "${SWIG_DEFS}")

include(swig_python_tools)

#swig_module_setup(${COMPONENT}_PYTHON_MODULES)

#build_python_tests()