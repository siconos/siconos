#================================================================
# cmake utilities to build kernel component
#================================================================
#
# The objective is to call component_setup to create the target <COMPONENT>.
# Before, it's necessary to set:
# 
# - COMPONENT component name
# - <COMPONENT>_DIRS: the list of paths (relative to CMAKE_CURRENT_SOURCE_DIR) that
#   contain source files
# - <COMPONENT>_EXCLUDE_SRCS: the list of files, in <COMPONENT>_DIRS, that must be excluded
#   from build.
# - <COMPONENT>_INTERFACE_INCLUDE_DIRECTORIES: a list of directories
#   to populate the interface of the target for include directories at build time

# Component name (i.e. target name)
set(COMPONENT kernel)
message("-- Set up for ${PROJECT_NAME}_${COMPONENT} library ...\n")

# ------ source directories for current component ------
# What is needed by component to compile ?
# List here all directories that contain sources files
# for current component.
# Path must be relative to component path (i.e. to CMAKE_CURRENT_SOURCE_DIR)
set(${COMPONENT}_DIRS
  src/.
  src/utils
  src/utils/SiconosAlgebra
  src/utils/SiconosMemory
  src/utils/SiconosSharedLibrary
  src/utils/SiconosTools
  src/utils/ProgressBar
  src/plugin
  src/modelingTools
  src/simulationTools
)

# -- Documentation --
# List of directories for which no doxygen doc will be generated
# By default all directories matching "test" are excluded.
set(${COMPONENT}_EXCLUDE_DOXY src/utils/SiconosTools
  # if included, it breaks swig process for mechanics
  )

# ------ include interface ------
# What is needed at build time
# by other targets to compile with current component.
# 
# It means that a call to
#  target_link_libraries(truc PRIVATE kernel)
# will imply -I<dirs> with dirs listed in
# ${COMPONENT}_INTERFACE_INCLUDE_DIRECTORIES.
#
# If all dirs are required, just set
# set(${COMPONENT}_INTERFACE_INCLUDE_DIRECTORIES ${${COMPONENT}_DIRS})
set(${COMPONENT}_INTERFACE_INCLUDE_DIRECTORIES ${${COMPONENT}_DIRS})

# ---- Final setup for the library ----

# Windows stuff
include(WindowsKernelSetup)

# -- create/setup component target --
include(ComponentSetup)
create_siconos_component(${COMPONENT})

# We require C++ 11. PUBLIC to propagate to other components.
if(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.8.0")
  target_compile_features(kernel PUBLIC cxx_std_11)
else()
  set_property(TARGET kernel PROPERTY CXX_STANDARD 11)
endif()

# --- Extra setup for the component ---

# Links with other Siconos components
target_link_libraries(kernel PRIVATE externals numerics)

# Links with non-Siconos libraries

# GMP
find_package(GMP REQUIRED)
target_link_libraries(kernel PRIVATE GMP::GMP)

# Boost must be set as a public dependency because of SiconosAlgebraTypeDefs.hpp
# This has to be reviewed !!!
target_link_libraries(kernel PUBLIC Boost::boost)

if(WITH_BOOST_LOG)
  find_package(Boost 1.61 REQUIRED COMPONENTS log)
  target_compile_definition(kernel PRIVATE BOOST_LOG_DYN_LINK)
endif()

# --- python bindings ---
if(WITH_${COMPONENT}_PYTHON_WRAPPER)
  add_subdirectory(swig)
endif()

if(WITH_PYTHON_WRAPPER)
  add_dependencies(${COMPONENT} ${COMPONENT}_docstrings)
endif()  

# ---- Installation ----
# Call siconos_component_install_setup(<COMPONENT>)
# to prepare installation of the current target.
#
# Before, it's necessary to set:
# 
# - <COMPONENT>_INSTALL_INTERFACE_INCLUDE_DIRECTORIES with all directories
#    that contain headers files that must be installed.
# 
set(${COMPONENT}_INSTALL_INTERFACE_INCLUDE_DIRECTORIES
  ${${COMPONENT}_DIRS} # All .hpp are installed
  )

# List header files (in dir above) that must not be installed 
#set(${COMPONENT}_HDRS_EXCLUDE)

siconos_component_install_setup(${COMPONENT})

# --- tests ---
include(${COMPONENT}_tests.cmake)


# CHECK THIS LATER :
# --- set linked libraries and linker language ---
#set(${COMPONENT}_LINKER_LANGUAGE CXX)
# CMAKE_DL_LIBS for dlopen or equivalent (used for plugins)
#list(APPEND ${COMPONENT}_LINK_LIBRARIES ${CMAKE_DL_LIBS})
