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
  src/bodies
  src/joints
  src/contactDetection/basicBroadphase
  )

if(WITH_OCC)
  list(APPEND ${COMPONENT}_DIRS src/occ)
endif()

if(WITH_BULLET)
  list(APPEND ${COMPONENT}_DIRS src/contactDetection/bullet)
endif()

if(WITH_MECHANISMS)
  list(APPEND ${COMPONENT}_DIRS
    mechanisms/CADMBTB
    mechanisms/MBTB
    mechanisms/TIMERMBTB
    mechanisms/frontEnd)
endif()

# Unstable sources. Will be included only
# if WITH_${COMPONENT}_UNSTABLE is true.
set(${COMPONENT}_Unstable_SRCS
  )

# List of directories of headers not to be installed
set(${COMPONENT}_HDRS_EXCLUDE_DIR)

# python_bindings - Set extra dependencies
if(WITH_${COMPONENT}_PYTHON_BINDINGS)
  set(SWIG_MODULE_${COMPONENT}_EXTRA_DEPS)
  # List of modules to be built in the python bindings
  set(${COMPONENT}_PYTHON_MODULES "Bodies;Joints")
  if(WITH_OCC)
    set(${COMPONENT}_PYTHON_MODULES "Bodies;Joints;OccWrap")
  endif()
endif()

include(WindowsMechanicsSetup)

# ---- Final setup for the library ----
# --- set linked libraries and linker language ---
set(${COMPONENT}_LINKER_LANGUAGE CXX)
set(${COMPONENT}_LINK_LIBRARIES ${SICONOS_LINK_LIBRARIES})
set(${COMPONENT}_LINK_LIBRARIES ${${COMPONENT}_LINK_LIBRARIES} kernel)
include(LibraryProjectSetup)
library_project_setup()

print_var(${COMPONENT}_SRCS)
# --- specific install ----


# if(WITH_${COMPONENT}_PYTHON_BINDINGS)
#   # Contact detection
#   install(PROGRAMS ${CMAKE_CURRENT_SOURCE_DIR}/swig_files/contact_detection/Tools.py
#   DESTINATION ${PYTHON_DIST_PACKAGES}/siconos/mechanics/contact_detection)

#   if(WITH_OCC)
#     install(PROGRAMS ${CMAKE_CURRENT_SOURCE_DIR}/swig_files/occ/__init__.py
#       DESTINATION ${SICONOS_PYTHON_INSTALL_DIR}}/mechanics/occ)
#   endif()
# endif()

# --- tests ---
include(${COMPONENT}_tests.cmake)
