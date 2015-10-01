#================================================================
# cmake utilities to build python-bindings for control component
#================================================================

# python_bindings - Set extra dependencies
set(SWIG_MODULE_${COMPONENT}_EXTRA_DEPS
  ${CMAKE_SOURCE_DIR}/Kernel/swig_files/Kernel.i
  ControlBase.i
  ControlTypemaps.i
  )

# List of modules to be built in the python bindings
set(${COMPONENT}_PYTHON_MODULES
  "${COMPONENT}/controller"
  "${COMPONENT}/observer"
  "${COMPONENT}/sensor"
  "${COMPONENT}/simulation"
)

set(${COMPONENT}_SWIG_INCLUDE_DIRECTORIES
  ${CMAKE_SOURCE_DIR}/Kernel/swig_files)
set(${COMPONENT}_SWIG_DEFS "${SWIG_DEFS}")

include(swig_python_tools)
swig_module_setup(${COMPONENT}_PYTHON_MODULES)
build_python_tests()