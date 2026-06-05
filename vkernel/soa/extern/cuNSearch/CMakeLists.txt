#cmake_minimum_required(VERSION 3.12)
#Requires cmake 3.12 for first class cuda support with visual studio

#project(cuNSearch LANGUAGES CXX CUDA)
# Visual studio solution directories.
enable_language(CUDA)
#set(CUDA_NVCC_EXECUTABLE /usr/local/cuda-13.0/bin/nvcc)
set_property(GLOBAL PROPERTY USE_FOLDERS on)


option(CUNSEARCH_USE_DOUBLE_PRECISION "Use double precision." OFF)


set(CUNSEARCH_USE_DOUBLE_PRECISION FALSE)

if(CUNSEARCH_USE_DOUBLE_PRECISION)
	message(STATUS "cuNSearch::Real = double")
else()
	message(STATUS "cuNSearch::Real = float")
endif(CUNSEARCH_USE_DOUBLE_PRECISION)

if(CUNSEARCH_USE_DOUBLE_PRECISION)
	add_compile_options(-DCUNSEARCH_USE_DOUBLE_PRECISION)
endif(CUNSEARCH_USE_DOUBLE_PRECISION)

OPTION(BUILD_AS_SHARED_LIBS "Build all the libraries as shared" ON)
if (BUILD_AS_SHARED_LIBS)
	if(WIN32)
		set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON CACHE BOOL "Export all symbols")
	endif(WIN32)	
endif (BUILD_AS_SHARED_LIBS)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
SET(CMAKE_DEBUG_POSTFIX "_d")

set (INCLUDE_HEADERS
	include/PointSet.h
	include/ActivationTable.h
	include/Common.h
	include/cuNSearch.h
)

set (HEADER_FILES
	src/Types.h
	src/cuNSearchDeviceData.h
	src/GridInfo.h
	src/NotImplementedException.h
	src/PointSetImplementation.h
	src/cuNSearchKernels.cuh
	src/helper_linearIndex.h
	src/helper_mortonCode.h
	Utils/cuda_helper.h
	Utils/Timing.h
	Utils/IDFactory.h
)

set (SOURCE_FILES
	src/PointSet.cu
	src/PointSetImplementation.cu
	src/cuNSearch.cu
	src/cuNSearchDeviceData.cu
	src/cuNSearchKernels.cu
	Utils/cuda_helper.cu
	Utils/Timing.cpp
	Utils/IDFactory.cpp
)

if(BUILD_AS_SHARED_LIBS)
	add_library(cuNSearch SHARED ${INCLUDE_HEADERS} ${HEADER_FILES} ${SOURCE_FILES})
else()
	add_library(cuNSearch STATIC ${INCLUDE_HEADERS} ${HEADER_FILES} ${SOURCE_FILES})
endif()

set_target_properties(cuNSearch PROPERTIES CUDA_ARCHITECTURES "${CMAKE_CUDA_ARCHITECTURES}")

include(GenerateExportHeader)
generate_export_header(cuNSearch
	BASE_NAME cuNSearch
	EXPORT_MACRO_NAME cuNSearch_EXPORT
	EXPORT_FILE_NAME ${CMAKE_CURRENT_BINARY_DIR}/cuNSearch/cuNSearch_export.h
	)

target_include_directories(cuNSearch PUBLIC
	"include"
	"Utils"
	${CUDA_INCLUDE_DIRS}
	${CMAKE_BINARY_DIR}/cuNSearch
)
target_link_libraries(cuNSearch PUBLIC ${CUDA_LIBRARIES})
target_compile_definitions(cuNSearch PUBLIC $<$<CONFIG:DEBUG>:DEBUG>)

list(APPEND INCLUDE_HEADERS
			${CMAKE_CURRENT_BINARY_DIR}/cuNSearch/cuNSearch_export.h)

install(FILES ${INCLUDE_HEADERS}
	DESTINATION include/)

install(TARGETS cuNSearch
	RUNTIME DESTINATION bin
	LIBRARY DESTINATION ${SICONOS_PYTHON_INSTALL_DIR}/nonos
	ARCHIVE DESTINATION ${SICONOS_PYTHON_INSTALL_DIR}/nonos
	)

option(BUILD_DEMO "Build example of how to use this library."
		ON)
if(BUILD_DEMO)
	add_subdirectory(demo)
endif(BUILD_DEMO)


unset(USE_DOUBLE_PRECISION CACHE)
