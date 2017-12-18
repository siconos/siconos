#================================================================
# cmake utilities to build numerics component
#================================================================

set(COMPONENT numerics)
message("-- Set up for ${PROJECT_NAME}_${COMPONENT} library ...\n")
# Name used for the target to be built
set(COMPONENT_LIBRARY_NAME ${PROJECT_NAME}_${COMPONENT})

# --- Numerics optional dependencies ---
include(gams_setup)


# --- Sources and headers ---
# source directories for numerics
set(${COMPONENT}_DIRS
  src/.
  src/tools
  src/tools/InterfaceToPathFerris
  src/tools/internal
  src/AVI
  src/LCP
  src/MLCP
  src/SOCP
  src/QP
  src/Relay
  src/FrictionContact
  src/FrictionContact/Generated
  src/NCP
  src/MCP
  src/GenericMechanical
  src/VI
)

# Static verifications (ACSL)
if(WITH_FRAMA_C)
  include(FramaC)
  file(GLOB frama_c_files src/FrictionContact/Generated/*.c)
  foreach(f ${frama_c_files})
    add_frama_c_test(${f} ENTRY main)
  endforeach()
  add_frama_c_test(${CMAKE_CURRENT_SOURCE_DIR}/src/tools/test/test_static_op3x3.c ENTRY main INCLUDES ${CMAKE_SOURCE_DIR}/externals/tools)
endif()

# Unstable sources. Will be included only
# if WITH_${COMPONENT}_UNSTABLE is true.
set(${COMPONENT}_Unstable_SRCS)

# List of directories of headers not to be installed
set(${COMPONENT}_HDRS_EXCLUDE_DIR)

# List of specific headers not to be installed
set(${COMPONENT}_HDRS_EXCLUDE
  src/FrictionContact/Generated/funcodegen.h
  )

# List of directories for which no doxygen doc will be generated
# By default all directories matching "test" are excluded.
set(${COMPONENT}_EXCLUDE_DOXY
  src/Unstable_or_deprecated)

# Windows stuff
include(WindowsNumericsSetup)

# ---- Final setup for the library ----
# --- set linked libraries and linker language ---
set(${COMPONENT}_LINKER_LANGUAGE C)
# but CXX may be needed with MlcpSimplex
if(MLCPSIMPLEX_FOUND) 
  set(${COMPONENT}_LINKER_LANGUAGE CXX)
endif()
set(${COMPONENT}_LINK_LIBRARIES externals ${SICONOS_LINK_LIBRARIES})
set(COMPONENT_TEST_LIB "numerics-test")
include(LibraryProjectSetup)
library_project_setup()

if(BUILD_AS_CPP)
  file(GLOB_RECURSE C_FILES ${CMAKE_CURRENT_SOURCE_DIR} *.c)
  set_source_files_properties(${C_FILES} PROPERTIES LANGUAGE CXX)
  set(${COMPONENT}_LINKER_LANGUAGE CXX)
endif()

# --- tests ---
include(${COMPONENT}_tests.cmake)
