#================================================================
# cmake utilities to build kernel component
#================================================================

set(COMPONENT kernel)
message("-- Set up for ${PROJECT_NAME}_${COMPONENT} library ...\n")
# Name used for the target to be built
set(COMPONENT_LIBRARY_NAME ${PROJECT_NAME}_${COMPONENT})

# --- Search component dependencies ---

# --- Sources and headers ---
# source directories for numerics
set(${COMPONENT}_DIRS
  src/.
  src/utils
  src/utils/numeric_bindings
  src/utils/contract
  src/utils/SiconosAlgebra
  src/utils/SiconosException
  src/utils/SiconosMemory
  src/utils/SiconosSharedLibrary
  src/utils/SiconosTools
  src/model
  src/plugin
  src/modelingTools
  src/simulationTools
)

# Unstable sources. Will be included only
# if WITH_${COMPONENT}_UNSTABLE is true.
set(${COMPONENT}_Unstable_SRCS
  src/simulationTools/GlobalFrictionContact.cpp
  src/simulationTools/MLCP2.cpp
  )

# List of directories of headers not to be installed
set(${COMPONENT}_HDRS_EXCLUDE_DIR src/SOL/lumod-c)

# List of directories for which no doxygen doc will be generated
# By default all directories matching "test" are excluded.
set(${COMPONENT}_EXCLUDE_DOXY
  src/utils/numeric_bindings src/utils/AtlasLocal
  )

# ---- Final setup for the library ----
# --- set linked libraries and linker language ---
set(${COMPONENT}_LINKER_LANGUAGE CXX)
set(${COMPONENT}_LINK_LIBRARIES ${SICONOS_LINK_LIBRARIES})
set(${COMPONENT}_LINK_LIBRARIES ${${COMPONENT}_LINK_LIBRARIES} numerics)

include(LibraryProjectSetup)
library_project_setup()

include(WindowsKernelSetup)

# --- specific install ----
install(FILES ${BOOST_HDRS} 
  DESTINATION include/${PROJECT_NAME}/boost/numeric/bindings/)


# --- tests ---
include(${COMPONENT}_tests.cmake)
