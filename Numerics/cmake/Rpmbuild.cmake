# http://mehdi.rabah.free.fr/Rpmbuild.cmake

# Rpmbuild : Create rpm packages for your projects and sub projects. Written by Mehdi Rabah.
# Was heavily inspired by UseDebian (Mathieu Malaterre) and UseRPMTools (TSP Team) modules

# need /usr/bin/rpmbuild

# Supported variables
# Summary: ${PACKAGE_DESCRIPTION_SUMMARY}
# Name: ${PACKAGE_NAME}
# Version: ${PACKAGE_VERSION}
# Release: ${PACKAGE_RELEASE}
# License: ${PACKAGE_LICENSE}
# Group: ${PACKAGE_GROUP}
# Source: ${PACKAGE_SOURCE}
# Packager: ${PACKAGE_MAINTAINER_NAME} ${PACKAGE_MAINTAINER_EMAIL}
#
# %description
# ${PACKAGE_DESCRIPTION}

FIND_PROGRAM(RPMBUILD
    NAMES rpmbuild
    PATHS "/usr/bin")

IF ( RPMBUILD )
    GET_FILENAME_COMPONENT(RPMBUILD_PATH ${RPMBUILD} ABSOLUTE)
    MESSAGE(STATUS "Found rpmbuild : ${RPMBUILD_PATH}")
    SET(RPMBUILD_FOUND "YES")
ELSE ( RPMBUILD ) 
    MESSAGE(STATUS "rpmbuild NOT found. rpm generation will not be available")
    SET(RPMBUILD_FOUND "NO")
ENDIF ( RPMBUILD )

# Main and only command of this module. 
MACRO(ADD_RPM RPM_NAME)

  SET(RPM_ROOTDIR ${CMAKE_BINARY_DIR}/RPM)
  
  FILE(MAKE_DIRECTORY ${RPM_ROOTDIR})
  FILE(MAKE_DIRECTORY ${RPM_ROOTDIR}/tmp)
  FILE(MAKE_DIRECTORY ${RPM_ROOTDIR}/BUILD)
  FILE(MAKE_DIRECTORY ${RPM_ROOTDIR}/RPMS)
  FILE(MAKE_DIRECTORY ${RPM_ROOTDIR}/SOURCES)
  FILE(MAKE_DIRECTORY ${RPM_ROOTDIR}/SPECS)
  FILE(MAKE_DIRECTORY ${RPM_ROOTDIR}/SRPMS)
  
  SET ( SPEC_FILE ${PROJECT_BINARY_DIR}/${PROJECT_NAME}.spec )
  
  
  # First choice for spec file : user defined variables 
  IF ("${ARGV1}" STREQUAL "")

      # Check if the mandatory variables are here
      # TODO
      IF( FALSE )
         message ( FATAL_ERROR "ADD_RPM command was not correctly configured for ${PROJECT_NAME}. See the documentation for more details" )
      ENDIF( FALSE )

      # Writing the spec file
      ADD_CUSTOM_COMMAND(
        OUTPUT ${SPEC_FILE} 
        COMMAND ${CMAKE_COMMAND} -E echo "Summary: ${PACKAGE_DESCRIPTION_SUMMARY}" > ${SPEC_FILE}
        COMMAND ${CMAKE_COMMAND} -E echo "Name: ${PACKAGE_NAME}" >> ${SPEC_FILE}
        COMMAND ${CMAKE_COMMAND} -E echo "Version: ${PACKAGE_VERSION}" >> ${SPEC_FILE}
        COMMAND ${CMAKE_COMMAND} -E echo "Release: ${PACKAGE_RELEASE}" >> ${SPEC_FILE}
        COMMAND ${CMAKE_COMMAND} -E echo "License: ${PACKAGE_LICENSE}" >> ${SPEC_FILE}
        COMMAND ${CMAKE_COMMAND} -E echo "Group: ${PACKAGE_GROUP}" >> ${SPEC_FILE}
        COMMAND ${CMAKE_COMMAND} -E echo "Source: ${PROJECT_NAME}" >> ${SPEC_FILE}
        COMMAND ${CMAKE_COMMAND} -E echo "Packager: ${PACKAGE_MAINTAINER_NAME}" " <${PACKAGE_MAINTAINER_EMAIL}>" >> ${SPEC_FILE}
        COMMAND ${CMAKE_COMMAND} -E echo "" >> ${SPEC_FILE}
        COMMAND ${CMAKE_COMMAND} -E echo "%description" >> ${SPEC_FILE}
        COMMAND ${CMAKE_COMMAND} -E echo "${PACKAGE_DESCRIPTION}" >> ${SPEC_FILE}
        COMMAND ${CMAKE_COMMAND} -E echo "" >> ${SPEC_FILE}

        COMMAND ${CMAKE_COMMAND} -E echo "%build" >> ${SPEC_FILE}
        COMMAND ${CMAKE_COMMAND} -E echo "cd ${PROJECT_BINARY_DIR}" >> ${SPEC_FILE}
        COMMAND ${CMAKE_COMMAND} -E echo "make" >> ${SPEC_FILE}
        COMMAND ${CMAKE_COMMAND} -E echo "" >> ${SPEC_FILE}
        
        COMMAND ${CMAKE_COMMAND} -E echo "%install" >> ${SPEC_FILE}
        COMMAND ${CMAKE_COMMAND} -E echo "cd ${PROJECT_BINARY_DIR}" >> ${SPEC_FILE}
        COMMAND ${CMAKE_COMMAND} -E echo "make install DESTDIR=${RPM_ROOTDIR}/tmp" >> ${SPEC_FILE}
        COMMAND ${CMAKE_COMMAND} -E echo "" >> ${SPEC_FILE}

        COMMAND ${CMAKE_COMMAND} -E echo "%files" >> ${SPEC_FILE}
        COMMAND ${CMAKE_COMMAND} -E echo "%defattr\(-,root,root,-\)" >> ${SPEC_FILE}
        COMMAND ${CMAKE_COMMAND} -E echo "%dir ${CMAKE_INSTALL_PREFIX}" >> ${SPEC_FILE}
        COMMAND ${CMAKE_COMMAND} -E echo "\"${CMAKE_INSTALL_PREFIX}/*\"" >> ${SPEC_FILE}
        COMMAND ${CMAKE_COMMAND} -E echo "" >> ${SPEC_FILE}

        COMMENT   "Generating spec file"
        VERBATIM
      )   

  ELSE ("${ARGV1}" STREQUAL "")
     ADD_CUSTOM_COMMAND(
        OUTPUT ${SPEC_FILE}
        COMMAND  ${CMAKE_COMMAND} -E copy "${ARGV1}" ${SPEC_FILE}
        COMMENT  "Copying user specified spec file : ${ARGV1}"
      )
  ENDIF("${ARGV1}" STREQUAL "")
  
  # the final target:
  ADD_CUSTOM_TARGET ( ${RPM_NAME}_rpm
    COMMAND   ${RPMBUILD_PATH} -bb --define=\"_topdir ${RPM_ROOTDIR}\" --buildroot=${RPM_ROOTDIR}/tmp ${SPEC_FILE}
    DEPENDS   ${SPEC_FILE}
    COMMENT   "Generating rpm binary package" )
    
  ADD_CUSTOM_COMMAND( TARGET ${RPM_NAME}_rpm PRE_BUILD
    COMMAND rm -f ${SPEC_FILE} ) 
    
  
ENDMACRO(ADD_RPM RPM_NAME)
