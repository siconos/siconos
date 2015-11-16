MACRO(windows_library_extra_setup LIB_NAME TARGET)

if(BUILD_SHARED_LIBS)
  if(MSVC)
    find_program(CMAKE_NM NAMES ${_CMAKE_TOOLCHAIN_PREFIX}nm HINTS ${_CMAKE_TOOLCHAIN_LOCATION})
    include(Platform/Windows-GNU) # for proper prefixes and suffixes
    GET_FILENAME_COMPONENT(SHARED_LIB_WE
      lib${LIB_NAME} NAME_WE)
    SET_TARGET_PROPERTIES(${TARGET} PROPERTIES
      ENABLE_EXPORT 1
      LINK_FLAGS /DEF:"${CMAKE_CURRENT_BINARY_DIR}/${SHARED_LIB_WE}.def")
    IF(CROSSCOMPILING_LINUX_TO_WINDOWS)
      ADD_CUSTOM_COMMAND(TARGET ${TARGET}
        PRE_BUILD
        COMMAND ${CMAKE_NM} ARGS @${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/${TARGET}.dir/objects1.rsp
        | egrep ' "(D|T|B)" ' | cut -f 3 -d ' ' | sed  's,^_,,' | sed '1iEXPORTS' > ${SHARED_LIB_WE}.def
        ) # gruik gruik
    ELSE()
      ADD_CUSTOM_COMMAND(TARGET ${TARGET}
        PRE_BUILD
        COMMAND ${CMAKE_NM} ARGS @${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/${TARGET}.dir/objects1.rsp
        | bash ${CMAKE_SOURCE_DIR}/cmake/export_filter.sh > ${SHARED_LIB_WE}.def
        ) # gruik gruik
    ENDIF()
  ENDIF(MSVC)
ENDIF(BUILD_SHARED_LIBS)

ENDMACRO()
