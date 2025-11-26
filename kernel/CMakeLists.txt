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
  src/utils/SiconosTools
  src/utils/ProgressBar
  src/utils/plugins
  src/modelingTools
  src/simulationTools
)

# -- Documentation --
# List of directories for which no doxygen doc will be generated
# By default all directories matching "test" are excluded.
set(${COMPONENT}_EXCLUDE_DOXY src/utils/SiconosTools
  # if included, it breaks swig process for mechanics
  )

# ---- Final setup for the library ----

# Windows stuff
include(WindowsKernelSetup)

# -- create/setup component target --
include(ComponentSetup)

create_siconos_component(${COMPONENT})
foreach(dir IN LISTS kernel_DIRS)
  target_include_directories(kernel PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/${dir}>)
endforeach()

# We require C++ 20. PUBLIC to propagate to other components.
target_compile_features(kernel PUBLIC cxx_std_20)

# --- Extra setup for the component ---
# remove boring warnings from boost bindings includes
#set_source_files_properties(
  #${CMAKE_SOURCE_DIR}/kernel/src/utils/SiconosAlgebra/EigenProblems.cpp
  # PROPERTIES COMPILE_OPTIONS "-Wno-unused-variable;-Wno-unused-local-typedefs")



# Links with other Siconos components
target_link_libraries(kernel PUBLIC externals numerics)
if(WITH_FORTRAN)
  target_link_libraries(kernel PRIVATE externals_fortran) # lsodar and hem5
endif()
# boost bindings
target_include_directories(kernel PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/externals/numeric_bindings>)
target_include_directories(kernel PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/externals/boost_contribs>)
# Mostly for the fortran API (empty if HAS_FORTRAN is off but required anyway)
target_include_directories(kernel PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/externals/tools>)

# Links with non-Siconos libraries

# GMP
find_package(GMP REQUIRED)
target_link_libraries(kernel PRIVATE GMP::GMP)

# Boost must be set as a public dependency because of SiconosAlgebraTypeDefs.hpp
# This has to be reviewed !!!
target_link_libraries(kernel PUBLIC Boost::boost)
if(WITH_BOOST_LOG)
  find_package(Boost 1.61 REQUIRED COMPONENTS log)
  target_compile_definitions(kernel PRIVATE BOOST_LOG_DYN_LINK)
endif()

if(WITH_JSON OR WITH_JSON_INSTALL)
 include(json_nlohmann_setup)
  # Do not use implicit conversion from a json file
  # This is a recommandation from nlohmann doc.
 target_compile_definitions(kernel PUBLIC JSON_ImplicitConversions=OFF)
 target_link_libraries(kernel PUBLIC nlohmann_json::nlohmann_json)
endif()

# Eigen
find_package(Eigen3 REQUIRED)
# Should we propagate Eigen dep? Check this later ...
target_link_libraries (kernel PUBLIC Eigen3::Eigen) 
target_compile_definitions(kernel PUBLIC EIGEN_DEFAULT_DENSE_INDEX_TYPE=int) # To be complient with scipy sparse
target_compile_definitions(kernel PUBLIC SICONOS_SPARSE_STORAGE_INDEX_TYPE=int)
# Sparse storage
if(SICONOS_SPARSE)
  message("Activate sparse storage (Eigen) for all matrices in Siconos. Warning: expermimental setup, still under test!")
  target_compile_definitions(kernel PUBLIC SICONOS_SPARSE)
endif()

get_target_property(raw_libs kernel INTERFACE_LINK_LIBRARIES)
message("KERNEL TEST ${raw_libs}")


# --- python bindings ---
if(WITH_${COMPONENT}_PYTHON_WRAPPER)
  add_subdirectory(swig)
endif()
if(WITH_PYB11_WRAPPER)
  add_subdirectory(python)
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

siconos_component_install_setup(${COMPONENT})

# --- tests ---
include(${COMPONENT}_tests.cmake)


# CHECK THIS LATER :
# --- set linked libraries and linker language ---
#set(${COMPONENT}_LINKER_LANGUAGE CXX)
# CMAKE_DL_LIBS for dlopen or equivalent (used for plugins)
#list(APPEND ${COMPONENT}_LINK_LIBRARIES ${CMAKE_DL_LIBS})
