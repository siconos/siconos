
set(installed_targets ${installed_targets}
  CACHE INTERNAL "List of installed libraries for the siconos project.")

# List of cmake macros that will be distributed with siconos software.
# They may be required during call to siconos script
# or to configure projects depending on Siconos (e.g. siconos-tutorials)
set(cmake_macros
  SiconosTools.cmake
  FindPythonModule.cmake
  valgrind.supp
  )

foreach(file IN LISTS cmake_macros)
  install(FILES cmake/${file} DESTINATION ${ConfigPackageLocation})
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
    COMMAND cat ${CMAKE_CURRENT_BINARY_DIR}/python_install_manifest.txt >> ${CMAKE_CURRENT_BINARY_DIR}/install_manifest.txt
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
  install(FILES ${CMAKE_BINARY_DIR}/scripts/CMakeLists.txt DESTINATION ${ConfigPackageLocation})
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
  INSTALL_DESTINATION ${ConfigPackageLocation})

install(FILES ${CMAKE_BINARY_DIR}/SiconosConfig.h DESTINATION include/${PROJECT_NAME})

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
  DESTINATION ${ConfigPackageLocation}) 

# install config files
install(
  FILES ${CMAKE_BINARY_DIR}/siconos-config.cmake ${CMAKE_BINARY_DIR}/siconos-config-version.cmake
  DESTINATION ${ConfigPackageLocation})
