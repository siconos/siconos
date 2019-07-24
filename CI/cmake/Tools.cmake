# Get a list of uninitialized cache variables and append
# them to vars.
macro(get_uninitialized_vars vars)
  get_cmake_property(CACHE_VARS CACHE_VARIABLES)
  foreach(CACHE_VAR ${CACHE_VARS})
    get_property(CACHE_VAR_TYPE CACHE ${CACHE_VAR} PROPERTY TYPE)
    if(CACHE_VAR_TYPE STREQUAL "UNINITIALIZED")
      list(APPEND ${vars} -D${CACHE_VAR}=${${CACHE_VAR}})
    endif()
  endforeach()
endmacro()

# set an option from script
macro(set_option var value)
  list(APPEND SICONOS_CMAKE_OPTIONS -D${var}=${value})
endmacro()

# set -DSICONOS_COMPONENTS=...
macro(set_components)
  set(SICONOS_COMPONENTS)
  message("START ....")
  foreach(comp  ${ARGN})
	set(SICONOS_COMPONENTS "${SICONOS_COMPONENTS}\;${comp}")
  endforeach()
endmacro()
