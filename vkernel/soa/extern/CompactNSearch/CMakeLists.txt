cmake_minimum_required(VERSION 3.10)

project(CompactNSearch)

# Visual studio solution directories.
set_property(GLOBAL PROPERTY USE_FOLDERS on)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
	if (CMAKE_VERSION VERSION_GREATER "3.8")
		link_libraries(OpenMP::OpenMP_CXX)
	else()
		set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	endif()
endif(OPENMP_FOUND)

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

target_include_directories(CompactNSearch PUBLIC include)

install(FILES "include/CompactNSearch" ${HEADER_FILES}
	DESTINATION include/)

install(TARGETS CompactNSearch
	RUNTIME DESTINATION bin
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib)
	
OPTION(USE_DOUBLE_PRECISION "Use double precision"	ON)
if (USE_DOUBLE_PRECISION)
	add_definitions( -DUSE_DOUBLE)	
endif (USE_DOUBLE_PRECISION)	

if(WIN32)
	add_definitions( -D_SILENCE_STDEXT_ARR_ITERS_DEPRECATION_WARNING)
endif() 

option(BUILD_DEMO "Build example of how to use this library."
		ON)
if(BUILD_DEMO)
	add_subdirectory(demo)
endif(BUILD_DEMO)

