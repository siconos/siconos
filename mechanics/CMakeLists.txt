#================================================================
# cmake utilities to build mechanics component
#================================================================

# This component has two optional sub-parts:
# - collision/bullet which requires bullet physics software  (and thus -DWITH_BULLET=ON cmake option)
# - occ which requires open-cascade community edition (and thus -DWITH_OCE=ON cmake option)
# Those parts are automatically built and installed if bullet and/or oce are available and properly detected by cmake.


set(COMPONENT mechanics)
message("-- Set up for ${PROJECT_NAME}_${COMPONENT} library ...\n")

# Name used for the target to be built
set(COMPONENT_LIBRARY_NAME ${PROJECT_NAME}_${COMPONENT})


# --- Sources and headers for 'native' mechanics part ---
set(${COMPONENT}_DIRS
  src/.
  src/collision
  src/collision/native
  src/collision/native/bodies
  src/joints
  )

# --- Optional parts ---

# -- Bullet --
# Must :
# - check if bullet is available on the system
# Will:
# - add collision/bullet to the build process
# - set SICONOS_HAS_BULLET var (distributed in siconosConfig.cmake)
if(WITH_BULLET)
  
  find_package(Bullet REQUIRED)
  # --> set BULLET_INCLUDE_DIRS, BULLET_FOUND, BULLET_LIBRARIES.
  # --> Accept : BULLET_ROOT (as install path or Windows build path)
  if(BULLET_FOUND)
    set(SICONOS_HAS_BULLET TRUE CACHE INTERNAL "True if Bullet API has been found and is activated.")
    # If a custom bullet was set, provide it for user programs.
    message(STATUS "Bullet found.")
    message("    Bullet libraries : ${BULLET_LIBRARIES}.")
    message("    Bullet headers path : ${BULLET_INCLUDE_DIRS}.")
    # Append mechanics bullet component sources to sources list
    list(APPEND ${COMPONENT}_DIRS src/collision/bullet)
    # # Add bullet headers to search path
    # # 'private' add limited to current target would be better.
    # # But those headers are missing during swig process. Todo : find a proper
    # # way to deal with this.
    # # target_include_directories(${COMPONENT} PRIVATE ${BULLET_INCLUDE_DIRS})
    include_directories(${BULLET_INCLUDE_DIRS})
    # # Add bullet to component list of dep
    set(${COMPONENT}_LINK_LIBRARIES ${${COMPONENT}_LINK_LIBRARIES} ${BULLET_LIBRARIES})
  endif()
endif()



# -- OCE --
# Must :
# - check if oce is available on the system
# Will:
# - add occ to the build process
# - set SICONOS_HAS_OCE var (distributed in siconosConfig.cmake)
if(WITH_OCE)
  include(oce_setup)
  # Source dir for 'oce' component
  list(APPEND ${COMPONENT}_DIRS src/occ)
  # Add oce headers to search path
  # 'private' add limited to current target would be better.
  # But those headers are missing during swig process. Todo : find a proper
  # way to deal with this.
  # target_include_directories(${COMPONENT} PRIVATE ${OCE_INCLUDE_DIRS})
  include_directories(${OCE_INCLUDE_DIRS})
  # Add oce libs to components links
  set(${COMPONENT}_LINK_LIBRARIES ${${COMPONENT}_LINK_LIBRARIES} ${OCE_LIBRARIES})

endif()

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
list(APPEND ${COMPONENT}_LINK_LIBRARIES numerics externals kernel)




if(WITH_OCE OR CMAKE_SKIP_RPATH)
  # if no RPATH, then linking does not work for tests without specifying externals
  # for n2qn1
  set(${COMPONENT}_LINK_LIBRARIES ${${COMPONENT}_LINK_LIBRARIES} externals)
endif()

include(LibraryProjectSetup)

library_project_setup()


if(SICONOS_HAS_BULLET)
  # add bullet dep to siconos-mechanics interface (headers, libs and so on)
  #target_link_libraries(${COMPONENT} PUBLIC ${BULLET_LIBRARIES})

  # Update compile flags for component
  if(BULLET_USE_DOUBLE_PRECISION)
    # Do we need this as public or private? Check if mechanics io needs this option?
    target_compile_options(${COMPONENT} PUBLIC "-DBT_USE_DOUBLE_PRECISION")
  endif()
endif()


# --- tests ---
include(${COMPONENT}_tests.cmake)
