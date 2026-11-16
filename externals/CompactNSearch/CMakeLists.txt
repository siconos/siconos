cmake_minimum_required(VERSION 3.10)

project(CompactNSearch)

# Visual studio solution directories.
set_property(GLOBAL PROPERTY USE_FOLDERS on)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# OpenMP is linked PRIVATE below (local change): the headers do not use it,
# and a public OpenMP::OpenMP_CXX gave -fopenmp to every consumer of
# Siconos::vkernel, kernel-only programs included.
find_package(OpenMP REQUIRED)

OPTION(BUILD_AS_SHARED_LIBS "Build all the libraries as shared" OFF)
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
target_link_libraries(CompactNSearch PRIVATE OpenMP::OpenMP_CXX)

# BUILD/INSTALL_INTERFACE: the target is exported (siconosTargets) so that
# Siconos::vkernel can link it; a bare relative path cannot be exported.
# The headers are installed flat into include/ (see the install rule below).
target_include_directories(CompactNSearch PUBLIC
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
	$<INSTALL_INTERFACE:include>)

install(FILES "include/CompactNSearch" ${HEADER_FILES}
	DESTINATION include/)

# EXPORT siconosTargets: the vkernel component links this library and exports
# it (Siconos::CompactNSearch), so that a consumer of Siconos::vkernel gets
# the broad phase without naming the library itself.  Local addition.
install(TARGETS CompactNSearch
	EXPORT siconosTargets
	RUNTIME DESTINATION bin
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib)

OPTION(COMPACTNSEARCH_USE_DOUBLE_PRECISION "Use double precision"	OFF)
set(COMPACTNSEARCH_USE_DOUBLE_PRECISION OFF)
if (COMPACTNSEARCH_USE_DOUBLE_PRECISION)
	add_definitions(CompactNSearch -DCOMPACTNSEARCH_USE_DOUBLE)
endif (COMPACTNSEARCH_USE_DOUBLE_PRECISION)

if(WIN32)
	add_definitions(CompactNSearch -D_SILENCE_STDEXT_ARR_ITERS_DEPRECATION_WARNING)
endif()

option(BUILD_DEMO "Build example of how to use this library."
		OFF)
if(BUILD_DEMO)
	add_subdirectory(demo)
endif(BUILD_DEMO)
