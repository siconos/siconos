
set(installed_targets ${installed_targets}
  CACHE INTERNAL "List of installed libraries for the siconos project.")

# List of cmake macros that will be distributed with siconos software.
# They may be required during call to siconos script
# or to configure projects depending on Siconos (e.g. siconos-tutorials)
set(cmake_macros
  SiconosTools.cmake
  FindPythonModule.cmake
  valgrind.supp
  FindBLASDEV.cmake
  FindLAPACKDEV.cmake
  BlasLapackUtils.cmake
)


if(SICONOS_HAS_OpenCASCADE)
  list(APPEND cmake_macros occ_setup.cmake)
endif()

foreach(file IN LISTS cmake_macros)
  install(FILES cmake/${file} DESTINATION ${SiconosConfigPackageLocation})
endforeach()

# =========== uninstall target ===========
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/cmake_uninstall.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
  IMMEDIATE @ONLY)

if(WITH_PYTHON_WRAPPER)
  # deal with files installed for python 
  add_custom_target(uninstall
    echo >> ${CMAKE_CURRENT_BINARY_DIR}/install_manifest.txt
    #COMMAND cat ${CMAKE_CURRENT_BINARY_DIR}/python_install_manifest.txt >> ${CMAKE_CURRENT_BINARY_DIR}/install_manifest.txt
    COMMAND PYTHONPATH=${SICONOS_PYTHON_INSTALL_DIR} ${Python_EXECUTABLE} -m pip uninstall siconos
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
else()
  add_custom_target(uninstall
    echo >> ${CMAKE_CURRENT_BINARY_DIR}/install_manifest.txt
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
endif()


# ==== Generate siconos scripts ====
# Generate the driver file required to execute siconos.
if(EXISTS ${CMAKE_SOURCE_DIR}/scripts/CMakeLists.txt.in)
  configure_file(scripts/CMakeLists.txt.in scripts/CMakeLists-temp.txt @ONLY)
  configure_file(${CMAKE_BINARY_DIR}/scripts/CMakeLists-temp.txt ${CMAKE_BINARY_DIR}/scripts/CMakeLists.txt @ONLY)
  install(FILES ${CMAKE_BINARY_DIR}/scripts/CMakeLists.txt DESTINATION ${SiconosConfigPackageLocation})
endif()

if(EXISTS ${CMAKE_SOURCE_DIR}/scripts/siconos.py.in)
  message("-- Generate siconos script ...")
  configure_file(scripts/siconos.py.in scripts/siconos @ONLY)
  install(PROGRAMS ${CMAKE_BINARY_DIR}/scripts/siconos DESTINATION bin)
endif()

# ===== Siconos Package configuration file ====
# https://cmake.org/cmake/help/latest/manual/cmake-packages.7.html#creating-packages
# 
include(CMakePackageConfigHelpers)

# Generate ${PROJECT_NAME}Config.cmake
configure_package_config_file(siconos-config.cmake.in ${CMAKE_BINARY_DIR}/siconos-config.cmake
  INSTALL_DESTINATION ${SiconosConfigPackageLocation})

# Generate siconos-config-version.cmake file.
write_basic_package_version_file(
  "${CMAKE_BINARY_DIR}/siconos-config-version.cmake"
  VERSION ${MAJOR_VERSION}.${MINOR_VERSION}.${PATCH_VERSION}
  COMPATIBILITY ExactVersion
  )

export(EXPORT siconosTargets
  FILE "${CMAKE_CURRENT_BINARY_DIR}/siconosTargets.cmake"
  NAMESPACE Siconos::
  )

install(EXPORT siconosTargets
  NAMESPACE Siconos::
  DESTINATION ${SiconosConfigPackageLocation}) 

# install config files
install(
  FILES ${CMAKE_BINARY_DIR}/siconos-config.cmake ${CMAKE_BINARY_DIR}/siconos-config-version.cmake
  DESTINATION ${SiconosConfigPackageLocation})

if(WITH_GIT)
  # Save and install a file which contain git references for the current source directory
  # (branch and short commit number), mostly used by continuous integration and cdash 
  # to tag cdash builds.
  # git reference name (branch, tag ...) 
  execute_process(COMMAND
    ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
    OUTPUT_VARIABLE COMMIT_REF_NAME
    OUTPUT_STRIP_TRAILING_WHITESPACE
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
  file(WRITE ${CMAKE_BINARY_DIR}/siconos-commit.txt "${COMMIT_REF_NAME}-${SOURCE_ABBREV_GIT_SHA1}")
  install(FILES ${CMAKE_BINARY_DIR}/siconos-commit.txt DESTINATION ${SiconosConfigPackageLocation})
endif()
