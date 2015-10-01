if(BUILD_SHARED_LIBS)
  if(MSVC)
    find_program(CMAKE_NM NAMES ${_CMAKE_TOOLCHAIN_PREFIX}nm HINTS ${_CMAKE_TOOLCHAIN_LOCATION})
    include(Platform/Windows-GNU) # for proper prefixes and suffixes
    GET_FILENAME_COMPONENT(${COMPONENT}_SHARED_LIB_WE
      lib${COMPONENT_LIBRARY_NAME} NAME_WE)
    SET_TARGET_PROPERTIES(${COMPONENT} PROPERTIES
      ENABLE_EXPORT 1
      LINK_FLAGS /DEF:"${CMAKE_CURRENT_BINARY_DIR}/${${COMPONENT}_SHARED_LIB_WE}.def")
    IF(CROSSCOMPILING_LINUX_TO_WINDOWS)
      ADD_CUSTOM_COMMAND(TARGET ${COMPONENT}
        PRE_BUILD
        COMMAND ${CMAKE_NM} ARGS @${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/${COMPONENT}_shared.dir/objects1.rsp
        | egrep ' "(D|T|B)" ' | cut -f 3 -d ' ' | sed  's,^_,,' | sed '1iEXPORTS' > ${${COMPONENT}_SHARED_LIB_WE}.def
        ) # gruik gruik
    ELSE()
      ADD_CUSTOM_COMMAND(TARGET ${COMPONENT}
        PRE_BUILD
        COMMAND ${CMAKE_NM} ARGS @${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/${COMPONENT}_shared.dir/objects1.rsp
        | bash ${CMAKE_SOURCE_DIR}/cmake/export_filter.sh > ${${COMPONENT}_SHARED_LIB_WE}.def
        ) # gruik gruik
    ENDIF()
  ENDIF(MSVC)
ENDIF(BUILD_SHARED_LIBS)

