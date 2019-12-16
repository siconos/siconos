# set an option from script
macro(set_option var value)
  list(APPEND SICONOS_CMAKE_OPTIONS -D${var}=${value})
endmacro()

# set -DSICONOS_COMPONENTS=...
macro(set_components)
  set(SICONOS_COMPONENTS)
  foreach(comp  ${ARGN})
    set(SICONOS_COMPONENTS "${SICONOS_COMPONENTS}\;${comp}")
  endforeach()
endmacro()
