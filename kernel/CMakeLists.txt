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
  src/utils/SiconosAlgebra
  src/utils/SiconosException
  src/utils/SiconosMemory
  src/utils/SiconosSharedLibrary
  src/utils/SiconosTools
  src/plugin
  src/modelingTools
  src/simulationTools
)

set(${COMPONENT}_HDRS_EXCLUDE
  src/utils/SiconosTools/CxxStd.hpp
  )

# Unstable sources. Will be included only
# if WITH_${COMPONENT}_UNSTABLE is true.
set(${COMPONENT}_Unstable_SRCS
  )

# ---- Final setup for the library ----
# --- set linked libraries and linker language ---
set(${COMPONENT}_LINKER_LANGUAGE CXX)
# CMAKE_DL_LIBS for dlopen or equivalent (used for plugins)
list(APPEND ${COMPONENT}_LINK_LIBRARIES ${CMAKE_DL_LIBS})
list(APPEND ${COMPONENT}_LINK_LIBRARIES ${SICONOS_LINK_LIBRARIES})
list(APPEND ${COMPONENT}_LINK_LIBRARIES externals numerics)

include(WindowsKernelSetup)
# Some directories to exclude from xml to swig process
set(${COMPONENT}_EXCLUDE_DOXY
  src/utils/SiconosTools # if included, it breaks swig process for mechanics
  )

include(LibraryProjectSetup)
library_project_setup()

# --- dependencies ---


# boost
# check https://cmake.org/cmake/help/latest/module/FindBoost.html?highlight=boost
find_package(Boost 1.61 REQUIRED
  COMPONENTS
  OPTIONAL_COMPONENTS log)


# Boost must be set as a public dependency because of SiconosAlgebraTypeDefs.hpp
# This has to be reviewed !!!
target_link_libraries(kernel PUBLIC Boost::boost)
if(WITH_BOOST_LOG)
  target_compile_definition(kernel PRIVATE BOOST_LOG_DYN_LINK)
endif()


# blas 
target_link_libraries(kernel PRIVATE BLAS::BLAS)

# gmp
find_package(GMP REQUIRED)
target_link_libraries(kernel PRIVATE GMP::GMP)

# gmp
find_package(GMP REQUIRED)
target_link_libraries(kernel PRIVATE GMP::GMP)
# target_link_libraries(_kernel PRIVATE GMP::GMP) Note FP : failed due to swig setup. Fix this later




# --- tests ---
include(${COMPONENT}_tests.cmake)
