#
# Common setup
#

# an encourage to out of source builds 
INCLUDE(OutOfSourcesBuild)

INCLUDE(SiconosVersion)

# misc tools
INCLUDE(SiconosTools)

MACRO(SICONOS_PROJECT 
    _PROJECT_NAME)
  
  # Build options
  # Static and shared libs : defaults
  OPTION(BUILD_SHARED_LIBS "Building of shared libraries" ON)
  OPTION(BUILD_STATIC_LIBS "Building of static libraries" OFF)
  OPTION(WITH_TESTS_COVERAGE "Code coverage setup" OFF)
  OPTION(WITH_GIT "Consider sources are under GIT" OFF)
  OPTION(WITH_DEFAULT_BUILD_TYPE "Use a default build type (Release)" ON)
  OPTION(WITH_DOCUMENTATION "Build doxygen documentation with 'make doc'" OFF)
  OPTION(WITH_TESTING "Enable 'make test' target" ON)
  OPTION(WITH_TIMERS "Enable timers" OFF)
  OPTION(WITH_MUMPS "Compilation with MUMPS solver" OFF)
  OPTION(WITH_FCLIB "link with fclib when this mode is enable. Default = off." OFF)
  OPTION(WITH_SYSTEM_INFO "Print some CMake variables. Default = off." OFF)
  OPTION(WITH_CPACK "Configuration for cpack. Default = on." ON)

  # get system architecture 
  # https://raw.github.com/petroules/solar-cmake/master/TargetArch.cmake
  INCLUDE(TargetArch)
  TARGET_ARCHITECTURE(SYSTEM_ARCHITECTURE)

  # some informations
  IF(WITH_SYSTEM_INFO)
    INCLUDE(CMakePrintSystemInformation)
    MESSAGE(STATUS "SYSTEM ARCHITECTURE: ${SYSTEM_ARCHITECTURE}")
  ENDIF(WITH_SYSTEM_INFO)

  # features summary
  INCLUDE(FeatureSummary)

  # Build type
  IF(WITH_DEFAULT_BUILD_TYPE)
    IF(NOT CMAKE_BUILD_TYPE)
      SET(CMAKE_BUILD_TYPE "Release")
    ENDIF(NOT CMAKE_BUILD_TYPE)
  ENDIF(WITH_DEFAULT_BUILD_TYPE)

  # Project version
  STRING(TOLOWER ${_PROJECT_NAME} _LPROJECT_NAME)
  
  SET(PROJECT_SHORT_NAME ${_PROJECT_NAME})
  
  SET(PROJECT_PACKAGE_NAME "siconos-${_LPROJECT_NAME}")

  # Install lib directory 32, 64 etc. on Fedora, Debian 
  # http://public.kitware.com/Bug/view.php?id=11964
  INCLUDE(GNUInstallDirs)
  ASSERT(CMAKE_INSTALL_LIBDIR)

  # PACKAGE PROJECT SETUP
  PROJECT(${PROJECT_PACKAGE_NAME})

  SET(VERSION "${MAJOR_VERSION}.${MINOR_VERSION}.${PATCH_VERSION}")  
  
  # get git update command
  IF(WITH_GIT)
    FIND_PACKAGE(Git)
    MESSAGE(STATUS "git executable : ${GIT_EXECUTABLE}")
    MESSAGE(STATUS "git command : ${GITCOMMAND}")
    MESSAGE(STATUS "git update options : ${GIT_UPDATE_OPTIONS}")

    SET(CTEST_GIT_COMMAND "${GIT_EXECUTABLE}" )
    SET(UPDATE_COMMAND "${GITCOMMAND}")
    SET(UPDATE_OPTIONS "${GIT_UPDATE_OPTIONS}")

    EXECUTE_PROCESS(COMMAND 
      ${GIT_EXECUTABLE} log -n 1 --pretty=format:'%h' 
      OUTPUT_VARIABLE SOURCE_ABBREV_GIT_SHA1
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

    EXECUTE_PROCESS(COMMAND 
      ${GIT_EXECUTABLE} log -n 1 --pretty=format:'%H' 
      OUTPUT_VARIABLE SOURCE_GIT_SHA1
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
  ENDIF()
  # Some macros needed to check compilers environment
  INCLUDE(CheckSymbolExists)
  INCLUDE(CheckFunctionExists)
  INCLUDE(CheckIncludeFileCXX)
  INCLUDE(CheckIncludeFile)
  INCLUDE(CheckStructHasMember)

  # Compilers environment
  IF(CMAKE_C_COMPILER)
    INCLUDE(CheckCCompilerFlag)
    CHECK_C_COMPILER_FLAG("-std=c99" C_HAVE_C99)
    CHECK_C_COMPILER_FLAG("-Wall" C_HAVE_WALL)
    CHECK_C_COMPILER_FLAG("-lm" C_HAVE_LINKER_M)
    CHECK_C_COMPILER_FLAG("-Wextra -Wno-unused-parameter" C_HAVE_WEXTRA)
    CHECK_C_COMPILER_FLAG("-static -static-libgcc" C_HAVE_STATIC_LINK)
  ENDIF(CMAKE_C_COMPILER)

  IF(CMAKE_CXX_COMPILER)
    INCLUDE(TestCXXAcceptsFlag)
    CHECK_CXX_ACCEPTS_FLAG("-static -static-libgcc -static-libstdc++" CXX_HAVE_STATIC_LINK)
    CHECK_CXX_ACCEPTS_FLAG("-Woverloaded-virtual" CXX_HAVE_WOVERLOADED_VIRTUAL)
    CHECK_CXX_ACCEPTS_FLAG("-diag-disable 654" CXX_HAVE_DIAG_DISABLE_654)
    CHECK_CXX_ACCEPTS_FLAG("-D__aligned__=ignored" CXX_HAVE_D__ALIGNED__IGNORED)
  ENDIF(CMAKE_CXX_COMPILER)

  # Get c compiler version (cf FindBoost.cmake version 2.8.7)
  IF(CMAKE_C_COMPILER)
    EXEC_PROGRAM(${CMAKE_C_COMPILER}
      ARGS -dumpversion
      OUTPUT_VARIABLE C_COMPILER_VERSION)
  ENDIF(CMAKE_C_COMPILER)

  # Get fortran compiler version (cf FindBoost.cmake version 2.8.7)
  IF(CMAKE_Fortran_COMPILER)
    EXEC_PROGRAM(${CMAKE_Fortran_COMPILER}
      ARGS -dumpversion
      OUTPUT_VARIABLE Fortran_COMPILER_VERSION)
  ENDIF(CMAKE_Fortran_COMPILER)

  # Get cxx compiler version (cf FindBoost.cmake version 2.8.7)
  IF(CMAKE_CXX_COMPILER)
    EXEC_PROGRAM(${CMAKE_CXX_COMPILER}
      ARGS ${CMAKE_CXX_COMPILER_ARG1} -dumpversion
      OUTPUT_VARIABLE CXX_COMPILER_VERSION)
  ENDIF(CMAKE_CXX_COMPILER)

  # check some headers
  IF(CMAKE_SYSTEM_NAME MATCHES Windows)
    CHECK_STRUCT_HAS_MEMBER("struct timespec" tv_sec pthread.h HAVE_TIME_H)
  ELSE(CMAKE_SYSTEM_NAME MATCHES Windows)
    CHECK_STRUCT_HAS_MEMBER("struct timespec" tv_sec time.h HAVE_TIME_H)
  ENDIF(CMAKE_SYSTEM_NAME MATCHES Windows)
  CHECK_FUNCTION_EXISTS(gettimeofday HAVE_SYSTIMES_H)
  IF(MSVC)
    SET(BUILD_AS_CPP TRUE)
  ENDIF(MSVC)

  # Link external lib statically. This comes handy when we want to distribute
  # Siconos on Mac or Windows
  # TODO complete this for other lib (libxml2, gmp, boost, ...)
  IF(CROSSCOMPILING_LINUX_TO_WINDOWS)
    OPTION(LINK_STATICALLY "Link external libraries statically (on if crosscompiling from linux to windows)" ON)
  ENDIF(CROSSCOMPILING_LINUX_TO_WINDOWS)
  IF(LINK_STATICALLY)
    SET(BLA_STATIC TRUE) # For blas/lapack
    SET(Boost_USE_STATIC_LIBS TRUE)
    # For the compiler
    IF(NOT MSVC)
      IF(NOT (C_HAVE_STATIC_LINK AND CXX_HAVE_STATIC_LINK))
        message(FATAL_ERROR "Your compiler has to support static linking flags (-static -static-libgcc -static-libstdc++ -static-libgfortran)")
      ELSE()
        APPEND_C_FLAGS("-static -static-libgcc")
        APPEND_CXX_FLAGS("-static -static-libgcc -static-libstdc++")
      ENDIF()
    ENDIF(NOT MSVC)
    APPEND_Fortran_FLAGS("-static -static-libgcc -static-libgfortran") # XXX No test :( -- xhub
  ENDIF(LINK_STATICALLY)

  IF(CMAKE_SYSTEM_NAME MATCHES Windows)
    if (NOT MINGW)
      set(CMAKE_FIND_LIBRARY_PREFIXES "lib" "" ${CMAKE_FIND_LIBRARY_PREFIXES})
    endif()
    SET(EXE_EXT ".exe")
  ELSE()
    SET(EXE_EXT)
  ENDIF()

  # Some http://pipol.inria.fr configurations

  # system configuration directory
  SET(PIPOL_RC_DIR ${CMAKE_SOURCE_DIR}/../Build/Pipol)

  # specific cmake command
  SET(PIPOL_CONFIGURE_COMMAND cmake -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX} -DWITH_TESTING=True  -DWITH_GIT=0 -DSOURCE_ABBREV_GIT_SHA1=${SOURCE_ABBREV_GIT_SHA1} ${CMAKE_SOURCE_DIR})

  INCLUDE(Pipol)
  
  # Tests+Dashboard configuration
  IF(WITH_TESTING)
    IF(NOT BUILDNAME)
      SET(BUILDNAME "${_PROJECT_NAME}")
    ENDIF()
    INCLUDE(SiconosCTest)
    INCLUDE(CTest)
    ENABLE_TESTING()  # Useless? done in CTest.cmake
  
  ENDIF(WITH_TESTING)
  # The library build stuff
  INCLUDE(LibraryProjectSetup)
  
  # Doxygen documentation
  IF(WITH_DOCUMENTATION)
    INCLUDE(SiconosDoc)
  ENDIF(WITH_DOCUMENTATION)

  # Top level install
  SET(CMAKE_INCLUDE_CURRENT_DIR ON)
  INSTALL(FILES AUTHORS ChangeLog COPYING INSTALL README 
    DESTINATION share/doc/siconos-${VERSION}/${_PROJECT_NAME})
  

  # To find XXXConfig.h
  include_directories(${CMAKE_BINARY_DIR}) 
  # Sources
  IF(IS_DIRECTORY ${CMAKE_SOURCE_DIR}/src)
    add_subdirectory(src)
  ENDIF(IS_DIRECTORY ${CMAKE_SOURCE_DIR}/src)

  # NumericsConfig.h/KernelConfig.h generation
#  IF(EXISTS ${CMAKE_SOURCE_DIR}/config.h.cmake)
#    CONFIGURE_FILE(config.h.cmake ${PROJECT_SHORT_NAME}Config.h)
#  ENDIF(EXISTS ${CMAKE_SOURCE_DIR}/config.h.cmake)

  # man files
  IF(IS_DIRECTORY ${CMAKE_SOURCE_DIR}/man)
    CONFIGURE_FILE(man/siconos.1.in man/siconos.1)
    INSTALL(FILES ${CMAKE_BINARY_DIR}/man/siconos.1 DESTINATION share/man/man1)
  ENDIF(IS_DIRECTORY ${CMAKE_SOURCE_DIR}/man)

  # scripts
  FILE(GLOB SCRIPT_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} */*.script.in)
  IF(SCRIPT_FILES)
    FOREACH(_FILE ${SCRIPT_FILES})
      GET_FILENAME_COMPONENT(_BFILE ${_FILE} NAME_WE)
      GET_FILENAME_COMPONENT(_PFILE ${_FILE} PATH)
      FILE(MAKE_DIRECTORY ${_PFILE})
      CONFIGURE_FILE(${_FILE} ${_PFILE}/${_BFILE} @ONLY)
      INSTALL(FILES ${CMAKE_BINARY_DIR}/${_PFILE}/${_BFILE} DESTINATION bin RENAME ${_BFILE}
        PERMISSIONS OWNER_READ GROUP_READ WORLD_READ 
                    OWNER_EXECUTE GROUP_EXECUTE WORLD_EXECUTE)
    ENDFOREACH(_FILE ${SCRIPT_FILES})
  ENDIF(SCRIPT_FILES)

  # cmakelists
  FILE(GLOB CMAKELISTS_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} */CMakeLists.txt.in)
  IF(CMAKELISTS_FILES)
    FOREACH(_FILE ${CMAKELISTS_FILES})
      GET_FILENAME_COMPONENT(_BFILE ${_FILE} NAME_WE)
      GET_FILENAME_COMPONENT(_PFILE ${_FILE} PATH)
      FILE(MAKE_DIRECTORY ${_PFILE})
      CONFIGURE_FILE(${_FILE} ${_PFILE}/${_BFILE}.txt @ONLY)
      INSTALL(FILES ${CMAKE_BINARY_DIR}/${_PFILE}/${_BFILE}.txt DESTINATION share/${PROJECT_PACKAGE_NAME} RENAME ${_BFILE}.txt)
    ENDFOREACH(_FILE ${CMAKELISTS_FILES})
  ENDIF(CMAKELISTS_FILES)


#  configure_file(Config.cmake.in
#    "${CMAKE_BINARY_DIR}/Siconos${PROJECT_SHORT_NAME}Config.cmake")
#  configure_file(ConfigVersion.cmake.in
#    "${CMAKE_BINARY_DIR}/Siconos${PROJECT_SHORT_NAME}ConfigVersion.cmake" @ONLY)
#  install(FILES
#    "${CMAKE_BINARY_DIR}/Siconos${PROJECT_SHORT_NAME}Config.cmake"
#    "${CMAKE_BINARY_DIR}/Siconos${_PROJECT_NAME}ConfigVersion.cmake"
#    DESTINATION share/${PROJECT_PACKAGE_NAME})


  # xml
  IF(IS_DIRECTORY ${CMAKE_SOURCE_DIR}/config/xmlschema)
    FILE(GLOB _SFILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} config/xmlschema/*.xsd)
    FOREACH(_F ${_SFILES})
      CONFIGURE_FILE(${_F} ${CMAKE_CURRENT_BINARY_DIR}/${_F} COPYONLY)
      INSTALL(FILES ${_F} DESTINATION share/${PROJECT_PACKAGE_NAME})
    ENDFOREACH(_F ${_SFILES})
  ENDIF(IS_DIRECTORY ${CMAKE_SOURCE_DIR}/config/xmlschema)
  
  # Packaging

  IF(WITH_CPACK)
    INCLUDE(InstallRequiredSystemLibraries)
  
    SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "${PACKAGE_DESCRIPTION_SUMMARY}")
    SET(CPACK_PACKAGE_DESCRIPTION "${PACKAGE_DESCRIPTION}")
    SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/README")
    
    # a package generation failure on mac ...
    IF(APPLE)
    ELSE(APPLE)
      SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/COPYING")
    ENDIF(APPLE)
    
    IF(PIPOL_IMAGE)
      SET(CPACK_PACKAGE_FILE_NAME "${PROJECT_PACKAGE_NAME}-${MAJOR_VERSION}.${MINOR_VERSION}.${PATCH_VERSION}-r${SOURCE_ABBREV_GIT_SHA1}-${CMAKE_BUILD_TYPE}-${PIPOL_IMAGE_NAME}")
    ENDIF(PIPOL_IMAGE)
    SET(CPACK_PACKAGE_NAME "${PROJECT_PACKAGE_NAME}")
  
    SET(CPACK_PACKAGE_VERSION_MAJOR "${MAJOR_VERSION}")
    SET(CPACK_PACKAGE_VERSION_MINOR "${MINOR_VERSION}")
    SET(CPACK_PACKAGE_VERSION_PATCH "${PATCH_VERSION}")
    SET(CPACK_PACKAGE_INSTALL_DIRECTORY "Siconos-${MAJOR_VERSION}.${MINOR_VERSION}.${PATCH_VERSION}")
    
    SET(CPACK_PACKAGE_CONTACT "Siconos development team <siconos-team@lists.gforge.inria.fr>")
    
    SET(CPACK_DEBIAN_PACKAGE_DEPENDS ${DEBIAN_PACKAGE_DEPENDS})
  
    INCLUDE(CPack)

  ENDIF(WITH_CPACK)

ENDMACRO(SICONOS_PROJECT)

MACRO(CLOSE_PROJECT)
  # NumericsConfig.h/KernelConfig.h generation
  IF(EXISTS ${CMAKE_SOURCE_DIR}/config.h.cmake)
    CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/config.h.cmake ${CMAKE_BINARY_DIR}/${PROJECT_SHORT_NAME}Config.h)
  ENDIF(EXISTS ${CMAKE_SOURCE_DIR}/config.h.cmake)
  
  FEATURE_SUMMARY(WHAT ALL)

ENDMACRO(CLOSE_PROJECT)

MACRO(WRITE_NOTES)
  IF(IS_DIRECTORY ${CMAKE_BINARY_DIR}/Testing)
    # a note file for the dashboard
    FILE(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/Testing/Notes)
    FILE(WRITE ${CMAKE_BINARY_DIR}/Testing/Notes/Build "git sha1 : ${SOURCE_GIT_SHA1}\n")
    FILE(APPEND ${CMAKE_BINARY_DIR}/Testing/Notes/Build "cmake version : ${CMAKE_VERSION}\n")
    # the default buildname
    FILE(APPEND ${CMAKE_BINARY_DIR}/Testing/Notes/Build "System name : ${CMAKE_SYSTEM_NAME}\n")
    FILE(APPEND ${CMAKE_BINARY_DIR}/Testing/Notes/Build "Processor   : ${CMAKE_SYSTEM_PROCESSOR}\n")
    FILE(APPEND ${CMAKE_BINARY_DIR}/Testing/Notes/Build "C compiler : ${CMAKE_C_COMPILER}\n")
    FILE(APPEND ${CMAKE_BINARY_DIR}/Testing/Notes/Build "C compiler version : ${C_COMPILER_VERSION}\n")
    FILE(APPEND ${CMAKE_BINARY_DIR}/Testing/Notes/Build "CXX compiler : ${CMAKE_CXX_COMPILER}\n")
    FILE(APPEND ${CMAKE_BINARY_DIR}/Testing/Notes/Build "CXX compiler version : ${CXX_COMPILER_VERSION}\n")
    FILE(APPEND ${CMAKE_BINARY_DIR}/Testing/Notes/Build "Fortran compiler : ${CMAKE_Fortran_COMPILER}\n")
    FILE(APPEND ${CMAKE_BINARY_DIR}/Testing/Notes/Build "Fortran compiler version : ${Fortran_COMPILER_VERSION}\n")
    FILE(APPEND ${CMAKE_BINARY_DIR}/Testing/Notes/Build "BLAS libraries : ${BLAS_LIBRARIES}\n")
    FILE(APPEND ${CMAKE_BINARY_DIR}/Testing/Notes/Build "LAPACK libraries : ${LAPACK_LIBRARIES}\n")
    FILE(APPEND ${CMAKE_BINARY_DIR}/Testing/Notes/Build "all libraries : ${${PROJECT_NAME}_LINK_LIBRARIES}\n")
  ENDIF(IS_DIRECTORY ${CMAKE_BINARY_DIR}/Testing)
ENDMACRO(WRITE_NOTES)
