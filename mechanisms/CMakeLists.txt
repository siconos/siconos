#================================================================
# cmake utilities to build mechanisms component
#================================================================

set(COMPONENT mechanisms)
message("-- Set up for ${PROJECT_NAME}_${COMPONENT} library ...\n")
# Name used for the target to be built
set(COMPONENT_LIBRARY_NAME ${PROJECT_NAME}_${COMPONENT})

# --- Search component dependencies ---
if(NOT WITH_OCE)
  message(FATAL_ERROR "Please run cmake again with option WITH_OCE=ON (required for mechanisms) or do not setup Siconos with mechanisms component.")
endif()

include(oce_setup)
# Add oce headers to search path
# 'private' add limited to current target would be better.
# But those headers are missing during swig process. Todo : find a proper
# way to deal with this.
# target_include_directories(${COMPONENT} PRIVATE ${OCE_INCLUDE_DIRS})
include_directories(${OCE_INCLUDE_DIRS})
# Add oce libs to components links
set(${COMPONENT}_LINK_LIBRARIES ${${COMPONENT}_LINK_LIBRARIES} ${OCE_LIBRARIES})

# --- Sources and headers ---
# source directories for numerics
set(${COMPONENT}_DIRS
  src/
  src/CADMBTB
  src/MBTB
  src/TIMERMBTB)

# Unstable sources. Will be included only
# if WITH_${COMPONENT}_UNSTABLE is true.
set(${COMPONENT}_Unstable_SRCS
  )

# List of directories of headers not to be installed
set(${COMPONENT}_HDRS_EXCLUDE_DIR)

# ---- Final setup for the library ----
# --- set linked libraries and linker language ---
set(${COMPONENT}_LINKER_LANGUAGE CXX)
list(APPEND ${COMPONENT}_LINK_LIBRARIES ${SICONOS_LINK_LIBRARIES})
list(APPEND ${COMPONENT}_LINK_LIBRARIES numerics externals kernel mechanics)


include(LibraryProjectSetup)
library_project_setup()

# --- tests ---
include(${COMPONENT}_tests.cmake)
