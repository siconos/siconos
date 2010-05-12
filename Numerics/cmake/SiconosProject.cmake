#
# Common setup
#

# before everything
CMAKE_MINIMUM_REQUIRED(VERSION 2.4.4)

# an encourage to out of source builds 
INCLUDE(OutOfSourcesBuild)

# misc tools
INCLUDE(SiconosTools)

MACRO(SICONOS_PROJECT 
    _PROJECT_NAME 
    MAJOR_VERSION MINOR_VERSION PATCH_VERSION)

  # Set cmake policies (cmake >= 2.6)
  IF(COMMAND CMAKE_POLICY)

    CMAKE_POLICY(VERSION 2.6.0)

    # minimum version required
    CMAKE_POLICY(SET CMP0000 NEW) 

    # CMAKE_BACKWARDS_COMPATIBILITY should no longer be used
    CMAKE_POLICY(SET CMP0001 NEW) 
    
    # logical target names must be globally unique
    CMAKE_POLICY(SET CMP0002 NEW) 
    
    # Libraries linked via full path no longer produce linker search
    # paths
    CMAKE_POLICY(SET CMP0003 NEW)

    # Libraries linked may not have leading or trailing white space
    CMAKE_POLICY(SET CMP0004 NEW) 
                               
    # Preprocessor definition values are now escaped automatically.
    CMAKE_POLICY(SET CMP0005 NEW)

    # Installing MACOSX_BUNDLE targets requires a BUNDLE DESTINATION.
    CMAKE_POLICY(SET CMP0006 NEW)
    
    #list command no longer ignores empty elements.
    CMAKE_POLICY(SET CMP0007 NEW)

  ENDIF(COMMAND CMAKE_POLICY)

  # Build options
  # Static and shared libs : defaults
  OPTION(BUILD_SHARED_LIBS "Building of shared libraries" ON)
  OPTION(BUILD_STATIC_LIBS "Building of static libraries" OFF)
  OPTION(WITH_TESTS_COVERAGE "Code coverage setup" OFF)
  OPTION(WITH_SVN "Consider SVN is online" OFF)
  OPTION(WITH_DEFAULT_BUILD_TYPE "Use a default build type (Release)" ON)
  OPTION(WITH_DOCUMENTATION "Build doxygen documentation with 'make doc'" OFF)
  OPTION(WITH_TESTING "Enable 'make test' target" ON)


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

  # Compilers environment
  IF(CMAKE_C_COMPILER)
    INCLUDE(CheckCCompilerFlag)
    CHECK_C_COMPILER_FLAG("-std=c99" C_HAVE_C99)
    CHECK_C_COMPILER_FLAG("-Wall" C_HAVE_WALL)
  ENDIF(CMAKE_C_COMPILER)

 # IF(CMAKE_CXX_COMPILER)
 #   INCLUDE(TestCXXAcceptsFlag)
#    CHECK_CXX_ACCEPTS_FLAG(-ffriend-injection CXX_HAVE_FRIEND_INJECTION)
 # ENDIF(CMAKE_CXX_COMPILER)

# IF(CMAKE_Fortran_COMPILER)
#    INCLUDE(fortran)
#    INCLUDE(FortranLibraries)
#  ENDIF(CMAKE_Fortran_COMPILER)

  # Some http://pipol.inria.fr configurations

  # system configuration directory
  SET(PIPOL_RC_DIR ${CMAKE_SOURCE_DIR}/../Build/Pipol)

  # specific cmake command
  SET(PIPOL_CONFIGURE_COMMAND cmake -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX} -DWITH_TESTING=True -DWITH_SVN=FALSE -DSVN_REVISION=${SVN_REVISION} ${CMAKE_SOURCE_DIR})

  INCLUDE(Pipol)

  # Tests+Dashboard configuration
  IF(WITH_TESTING)
    ENABLE_TESTING()
    
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
    
    INCLUDE(DartConfig)
    INCLUDE(Dart)

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
    SUBDIRS(src)
  ENDIF(IS_DIRECTORY ${CMAKE_SOURCE_DIR}/src)

  # To save build settings
  INCLUDE(CMakeExportBuildSettings)

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
  SET(CPACK_PACKAGE_INSTALL_DIRECTORY "CMake_${CMake_VERSION_MAJOR}.${CMake_VERSION_MINOR}")
  
  SET(CPACK_PACKAGE_CONTACT "Siconos development team")


  IF(PIPOL_IMAGE MATCHES "etch")
    SET(DEBIAN_PACKAGE_DEPENDS 
      "g++, gfortran, atlas3-base, atlas3-base-dev, atlas3-headers, libxml2-dev, libboost-dev, libboost-graph-dev, libgfortran1-dev, libgmp3-dev, siconos-numerics")
  ENDIF(PIPOL_IMAGE MATCHES "etch")

  IF(PIPOL_IMAGE MATCHES "lenny\\|ubuntu")
    SET(DEBIAN_PACKAGE_DEPENDS
      "g++, gfortran, libatlas-base-dev, libatlas-headers, libboost-dev, libboost-graph-dev, libgmp3-dev, libxml2-dev, siconos-numerics")
  ENDIF(PIPOL_IMAGE MATCHES "lenny\\|ubuntu")  

  SET(CPACK_DEBIAN_PACKAGE_DEPENDS ${DEBIAN_PACKAGE_DEPENDS})
  
  INCLUDE(CPack)

ENDMACRO(SICONOS_PROJECT)
