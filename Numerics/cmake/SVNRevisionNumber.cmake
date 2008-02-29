#
# SVN revision number
#
SET( ENV{LC_ALL} C )
IF(EXISTS ${CMAKE_SOURCE_DIR}/.svn)
  FIND_PACKAGE(Subversion)
  IF(Subversion_FOUND)
    Subversion_WC_INFO(${PROJECT_SOURCE_DIR} Project)
    IF(Project_WC_REVISION)
      SET(SVN_REVISION ${Project_WC_REVISION})
      MESSAGE(STATUS "This is a build from sources under svn")
      MESSAGE(STATUS "Current svn revision is ${SVN_REVISION}")
    ENDIF(Project_WC_REVISION)
  ENDIF(Subversion_FOUND)
ELSE(EXISTS ${CMAKE_SOURCE_DIR}/.svn)
  MESSAGE(STATUS "This is a build out of svn")
ENDIF(EXISTS ${CMAKE_SOURCE_DIR}/.svn)
