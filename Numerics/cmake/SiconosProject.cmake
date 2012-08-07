#
# Common setup
#

# an encourage to out of source builds 
INCLUDE(OutOfSourcesBuild)

# misc tools
INCLUDE(SiconosTools)

MACRO(SICONOS_PROJECT 
    _PROJECT_NAME 
    MAJOR_VERSION MINOR_VERSION PATCH_VERSION)
  
  # Build options
  # Static and shared libs : defaults
  OPTION(BUILD_SHARED_LIBS "Building of shared libraries" ON)
  OPTION(BUILD_STATIC_LIBS "Building of static libraries" OFF)
  OPTION(WITH_TESTS_COVERAGE "Code coverage setup" OFF)
  OPTION(WITH_SVN "Consider SVN is online" OFF)
  OPTION(WITH_DEFAULT_BUILD_TYPE "Use a default build type (Release)" ON)
  OPTION(WITH_DOCUMENTATION "Build doxygen documentation with 'make doc'" OFF)
  OPTION(WITH_TESTING "Enable 'make test' target" ON)
  OPTION(WITH_TIMERS "Enable timers" OFF)
  OPTION(WITH_MUMPS "Compilation with MUMPS solver" OFF)
  OPTION(WITH_FCLIB "link with fclib when this mode is enable. Default = off." OFF)

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
  
  # PACKAGE PROJECT SETUP
  PROJECT(${PROJECT_PACKAGE_NAME})

  SET(VERSION "${MAJOR_VERSION}.${MINOR_VERSION}.${PATCH_VERSION}")  
  
  # try to get the SVN revision number
  IF(WITH_SVN)
    INCLUDE(SVNRevisionNumber)
  ENDIF(WITH_SVN)

  # Some macros needed to check compilers environment
  INCLUDE(CheckSymbolExists)
  INCLUDE(CheckFunctionExists)
  INCLUDE(CheckIncludeFileCXX)
  INCLUDE(CheckIncludeFile)

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
  ENDIF(CMAKE_CXX_COMPILER)

  # check some headers
  CHECK_INCLUDE_FILE(time.h HAVE_TIME_H)
  CHECK_INCLUDE_FILE(sys/times.h HAVE_SYSTIMES_H)

  # Link external lib statically. This icomes handy when we want to distribute
  # Siconos on Mac or Windows
  # TODO complete this for other lib (libxml2, gmp, boost, ...)
  IF(CROSSCOMPILING_LINUX_TO_WINDOWS)
    OPTION(LINK_STATICALLY "Link external libraries statically (on if crosscompiling from linux to windows)" ON)
  ENDIF()
  IF(LINK_STATICALLY)
    SET(BLA_STATIC TRUE) # For blas/lapack
    # For the compiler
    IF(NOT (C_HAVE_STATIC_LINK AND CXX_HAVE_STATIC_LINK))
      message(FATAL_ERROR "Your compiler has to support static linking flags (-static -static-libgcc -static-libstdc++ -static-libgfortran)")
    ELSE()
      APPEND_C_FLAGS("-static -static-libgcc")
      APPEND_CXX_FLAGS("-static -static-libgcc -static-libstdc++")
      APPEND_Fortran_FLAGS("-static -static-libgcc -static-libgfortran") # XXX No test :( -- xhub
    ENDIF()
  ENDIF(LINK_STATICALLY)

  IF(CMAKE_SYSTEM_NAME MATCHES Windows)
    SET(EXE_EXT ".exe")
  ELSE()
    SET(EXE_EXT)
  ENDIF()

  # Some http://pipol.inria.fr configurations

  # system configuration directory
  SET(PIPOL_RC_DIR ${CMAKE_SOURCE_DIR}/../Build/Pipol)

  # specific cmake command
  SET(PIPOL_CONFIGURE_COMMAND cmake -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX} -DWITH_TESTING=True -DWITH_SVN=FALSE -DSVN_REVISION=${SVN_REVISION} ${CMAKE_SOURCE_DIR})

  INCLUDE(Pipol)

  # Tests+Dashboard configuration
  IF(WITH_TESTING)
    
    IF(IS_DIRECTORY ${CMAKE_BINARY_DIR}/Testing)
      # a note file for the dashboard
      FILE(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/Testing/Notes)
      FILE(WRITE ${CMAKE_BINARY_DIR}/Testing/Notes/Build "svn revision : ${SVN_REVISION}\n")
      # the default buildname
      FILE(APPEND ${CMAKE_BINARY_DIR}/Testing/Notes/Build "System name : ${CMAKE_SYSTEM_NAME}\n")
      FILE(APPEND ${CMAKE_BINARY_DIR}/Testing/Notes/Build "Processor   : ${CMAKE_SYSTEM_PROCESSOR}\n")
      FILE(APPEND ${CMAKE_BINARY_DIR}/Testing/Notes/Build "C compiler : ${CMAKE_C_COMPILER}\n")
      FILE(APPEND ${CMAKE_BINARY_DIR}/Testing/Notes/Build "CXX compiler : ${CMAKE_CXX_COMPILER}\n")
      FILE(APPEND ${CMAKE_BINARY_DIR}/Testing/Notes/Build "Fortran compiler : ${CMAKE_Fortran_COMPILER}\n")
    ENDIF(IS_DIRECTORY ${CMAKE_BINARY_DIR}/Testing)

    IF(BUILDNAME_OPTIONS)
      SET(BUILDNAME "${_PROJECT_NAME}-${BUILDNAME_OPTIONS}")
    ELSE(BUILDNAME_OPTIONS)
      SET(BUILDNAME "${_PROJECT_NAME}")
    ENDIF(BUILDNAME_OPTIONS)
    
    IF(CMAKE_BUILD_TYPE)
      SET(BUILDNAME "${BUILDNAME}-${CMAKE_BUILD_TYPE}")
    ENDIF(CMAKE_BUILD_TYPE)
    
    IF(PIPOL_IMAGE)
      SET(BUILDNAME "${BUILDNAME}-${PIPOL_IMAGE_NAME}")
      SET(SITE ${PIPOL_SITE})
    ENDIF(PIPOL_IMAGE)
    
    # Tests coverage (taken from ViSp)
    
    #
    # Note: all of this is done with a recent cmake version (>2.6.0) with:
    # cmake -DCMAKE_BUILD_TYPE=Profile
    #
    IF(WITH_TESTS_COVERAGE)
      # Add build options for test coverage. Currently coverage is only supported
      # on gcc compiler
      # Because using -fprofile-arcs with shared lib can cause problems like:
      # hidden symbol `__bb_init_func', we add this option only for static
      # library build
      SET(BUILD_SHARED_LIBS)
      SET(CMAKE_BUILD_TYPE Debug)
      CHECK_CXX_ACCEPTS_FLAG(-ftest-coverage CXX_HAVE_FTEST_COVERAGE)
      CHECK_CXX_ACCEPTS_FLAG(-fprofile-arcs CXX_HAVE_PROFILE_ARCS)
      CHECK_C_COMPILER_FLAG(-ftest-coverage C_HAVE_FTEST_COVERAGE)
      CHECK_C_COMPILER_FLAG(-fprofile-arcs C_HAVE_PROFILE_ARCS)
      IF(CXX_HAVE_FTEST_COVERAGE AND CXX_HAVE_PROFILE_ARCS)
        MESSAGE("Adding test coverage flags to CXX compiler : -ftest-coverage -fprofile-arcs")
        SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -ftest-coverage -fprofile-arcs")
        SET (CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -fprofile-arcs -ftest-coverage")
        SET (CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -fprofile-arcs -ftest-coverage")
      ENDIF(CXX_HAVE_FTEST_COVERAGE AND CXX_HAVE_PROFILE_ARCS)
      
      IF(C_HAVE_FTEST_COVERAGE)
        MESSAGE("Adding test coverage flags to C compiler : -ftest-coverage")
        SET(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -ftest-coverage -fprofile-arcs")
        SET (CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -fprofile-arcs -ftest-coverage")
        SET (CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -fprofile-arcs -ftest-coverage")
      ENDIF(C_HAVE_FTEST_COVERAGE)
      
    ENDIF(WITH_TESTS_COVERAGE)

    INCLUDE(CTest)
    ENABLE_TESTING()  # Useless? done in CTest.cmake
  
  ENDIF(WITH_TESTING)
  # The library build stuff
  INCLUDE(LibraryProjectSetup)
  
  # Doxygen documentation
  IF(WITH_DOCUMENTATION)
    INCLUDE(SiconosDoc)
  ENDIF(WITH_DOCUMENTATION)

  # NumericsConfig.h/KernelConfig.h and include
  IF(EXISTS ${CMAKE_SOURCE_DIR}/config.h.cmake)
    IF(NOT CONFIG_H_GLOBAL_CONFIGURED)
      SET(CONFIG_H_GLOBAL_CONFIGURED 1 CACHE BOOL "${PROJECT_SHORT_NAME}Config.h global generation." )
      CONFIGURE_FILE(config.h.cmake ${PROJECT_SHORT_NAME}Config.h)
    ENDIF(NOT CONFIG_H_GLOBAL_CONFIGURED)
    INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR})
  ENDIF(EXISTS ${CMAKE_SOURCE_DIR}/config.h.cmake)

  # Top level install
  SET(CMAKE_INCLUDE_CURRENT_DIR ON)
  INSTALL(FILES AUTHORS ChangeLog COPYING INSTALL README 
    DESTINATION share/doc/siconos-${VERSION}/${_PROJECT_NAME})

  # man files
  IF(IS_DIRECTORY ${CMAKE_SOURCE_DIR}/man)
    CONFIGURE_FILE(man/siconos.1.in man/siconos.1)
    INSTALL(FILES ${CMAKE_BINARY_DIR}/man/siconos.1 DESTINATION man/man1)
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


  # xml
  IF(IS_DIRECTORY ${CMAKE_SOURCE_DIR}/config/xmlschema)
    FILE(GLOB _SFILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} config/xmlschema/*.xsd)
    FOREACH(_F ${_SFILES})
      CONFIGURE_FILE(${_F} ${CMAKE_CURRENT_BINARY_DIR}/${_F} COPYONLY)
      INSTALL(FILES ${_F} DESTINATION share/${PROJECT_PACKAGE_NAME})
    ENDFOREACH(_F ${_SFILES})
  ENDIF(IS_DIRECTORY ${CMAKE_SOURCE_DIR}/config/xmlschema)
  
  # Sources
  IF(IS_DIRECTORY ${CMAKE_SOURCE_DIR}/src)
    add_subdirectory(src)
  ENDIF(IS_DIRECTORY ${CMAKE_SOURCE_DIR}/src)
  
  # Packaging
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
    SET(CPACK_PACKAGE_FILE_NAME "${PROJECT_PACKAGE_NAME}-${MAJOR_VERSION}.${MINOR_VERSION}.${PATCH_VERSION}-r${SVN_REVISION}-${CMAKE_BUILD_TYPE}-${PIPOL_IMAGE_NAME}")
  ENDIF(PIPOL_IMAGE)
  SET(CPACK_PACKAGE_NAME "${PROJECT_PACKAGE_NAME}")
  
  SET(CPACK_PACKAGE_VERSION_MAJOR "${MAJOR_VERSION}")
  SET(CPACK_PACKAGE_VERSION_MINOR "${MINOR_VERSION}")
  SET(CPACK_PACKAGE_VERSION_PATCH "${PATCH_VERSION}")
  SET(CPACK_PACKAGE_INSTALL_DIRECTORY "Siconos-${MAJOR_VERSION}.${MINOR_VERSION}.${PATCH_VERSION}")
  
  SET(CPACK_PACKAGE_CONTACT "Siconos development team <siconos-team@lists.gforge.inria.fr>")

  SET(CPACK_DEBIAN_PACKAGE_DEPENDS ${DEBIAN_PACKAGE_DEPENDS})
  
  INCLUDE(CPack)

ENDMACRO(SICONOS_PROJECT)
