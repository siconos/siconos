#================================================================
# cmake utilities to build mechanics component
#================================================================

set(COMPONENT mechanics)
message("-- Set up for ${PROJECT_NAME}_${COMPONENT} library ...\n")
# Name used for the target to be built
set(COMPONENT_LIBRARY_NAME ${PROJECT_NAME}_${COMPONENT})

# --- Search component dependencies ---

# --- Sources and headers ---
# source directories for numerics
set(${COMPONENT}_DIRS
  src/.
  src/collision
  src/collision/native
  src/collision/native/bodies
  src/joints
  )

if(WITH_OCC)
  list(APPEND ${COMPONENT}_DIRS src/occ)
endif()

if(WITH_MECHANISMS)
  list(APPEND ${COMPONENT}_DIRS
    src/mechanisms/CADMBTB
    src/mechanisms/MBTB
    src/mechanisms/TIMERMBTB)
endif()

# Unstable sources. Will be included only
# if WITH_${COMPONENT}_UNSTABLE is true.
set(${COMPONENT}_Unstable_SRCS
  )

# List of directories of headers not to be installed
set(${COMPONENT}_HDRS_EXCLUDE_DIR)

include(WindowsMechanicsSetup)

# ---- Final setup for the library ----
# --- set linked libraries and linker language ---
set(${COMPONENT}_LINKER_LANGUAGE CXX)
list(APPEND ${COMPONENT}_LINK_LIBRARIES ${SICONOS_LINK_LIBRARIES})
list(APPEND ${COMPONENT}_LINK_LIBRARIES numerics externals kernel)

if(WITH_BULLET)
  list(APPEND ${COMPONENT}_DIRS src/collision/bullet)
  include_directories(${BULLET_INCLUDE_DIRS})
  set(${COMPONENT}_LINK_LIBRARIES ${${COMPONENT}_LINK_LIBRARIES} ${BULLET_LIBRARIES})
endif()

if(WITH_MECHANISMS OR WITH_OCC OR CMAKE_SKIP_RPATH)
  # if no RPATH, then linking does not work for tests without specifying externals
  # for n2qn1
  set(${COMPONENT}_LINK_LIBRARIES ${${COMPONENT}_LINK_LIBRARIES} externals)
endif()

include(LibraryProjectSetup)
library_project_setup()

# --- tests ---
include(${COMPONENT}_tests.cmake)
