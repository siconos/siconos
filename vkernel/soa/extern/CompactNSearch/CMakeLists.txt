#cmake_minimum_required(VERSION 3.2)

#project(CompactNSearch)

# Visual studio solution directories.
set_property(GLOBAL PROPERTY USE_FOLDERS on)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if (UNIX)
find_package(OpenMP)
if (OPENMP_FOUND)
	if (CMAKE_VERSION VERSION_GREATER "3.8")
		link_libraries(OpenMP::OpenMP_CXX)
	else()
		set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	endif()
endif(OPENMP_FOUND)
endif (UNIX)

OPTION(BUILD_AS_SHARED_LIBS "Build all the libraries as shared" ON)
if (BUILD_AS_SHARED_LIBS)
	if(WIN32)
		add_definitions(-DCMAKE_WINDOWS_EXPORT_ALL_SYMBOLS)
		set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS YES CACHE BOOL "Export all symbols")
	endif(WIN32)
endif (BUILD_AS_SHARED_LIBS)

SET(CMAKE_DEBUG_POSTFIX "_d")

set (HEADER_FILES
		include/Config.h
		include/CompactNSearch.h
		include/PointSet.h
		include/DataStructures.h)

if (BUILD_AS_SHARED_LIBS)
	add_library(CompactNSearch SHARED
	    ${HEADER_FILES}
		src/CompactNSearch.cpp
)
else()
	add_library(CompactNSearch
	    ${HEADER_FILES}
		src/CompactNSearch.cpp
)
endif ()

target_include_directories(CompactNSearch PUBLIC include)

if (APPLE AND CMAKE_CXX_COMPILER_ID MATCHES "Clang")
	find_package(oneDPL REQUIRED)
	target_link_libraries(CompactNSearch PUBLIC oneDPL)
endif()

set_target_properties(CompactNSearch PROPERTIES
  INSTALL_RPATH "${SICONOS_PYTHON_INSTALL_DIR}/nonos"
  BUILD_WITH_INSTALL_RPATH TRUE
  INSTALL_RPATH_USE_LINK_PATH TRUE
)

install(FILES "include/CompactNSearch" ${HEADER_FILES}
	DESTINATION include/)

install(TARGETS CompactNSearch
	RUNTIME DESTINATION bin
	LIBRARY DESTINATION ${SICONOS_PYTHON_INSTALL_DIR}/nonos
	ARCHIVE DESTINATION ${SICONOS_PYTHON_INSTALL_DIR}/nonos)

OPTION(USE_DOUBLE_PRECISION "Use double precision"	ON)
if (USE_DOUBLE_PRECISION)
	add_definitions( -DUSE_DOUBLE)
endif (USE_DOUBLE_PRECISION)

option(BUILD_DEMO "Build example of how to use this library."
		ON)
if(BUILD_DEMO)
	add_subdirectory(demo)
endif(BUILD_DEMO)
