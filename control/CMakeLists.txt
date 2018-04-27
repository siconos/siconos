#================================================================
# cmake utilities to build control component
#================================================================

set(COMPONENT control)
message("-- Set up for ${PROJECT_NAME}_${COMPONENT} library ...\n")
# Name used for the target to be built
set(COMPONENT_LIBRARY_NAME ${PROJECT_NAME}_${COMPONENT})

# --- Search component dependencies ---

# --- Sources and headers ---
# source directories for numerics
set(${COMPONENT}_DIRS
  src/.
  src/utils
  src/Controller
  src/Observer
  src/Sensor
  src/Simulation
)

# Unstable sources. Will be included only
# if WITH_${COMPONENT}_UNSTABLE is true.
set(${COMPONENT}_Unstable_SRCS)

# List of directories of headers not to be installed
set(${COMPONENT}_HDRS_EXCLUDE_DIR)

# ---- Final setup for the library ----
# --- set linked libraries and linker language ---
set(${COMPONENT}_LINKER_LANGUAGE CXX)
list(APPEND ${COMPONENT}_LINK_LIBRARIES ${SICONOS_LINK_LIBRARIES})
list(APPEND ${COMPONENT}_LINK_LIBRARIES numerics kernel)
if (CMAKE_SKIP_RPATH)
  # if no RPATH, then linking does not work for tests without specifying externals
  list(APPEND ${COMPONENT}_LINK_LIBRARIES externals)
endif()

include(WindowsControlSetup)

include(LibraryProjectSetup)
library_project_setup()

# --- tests ---
include(${COMPONENT}_tests.cmake)
