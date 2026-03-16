file(GLOB GPU_SOURCES CONFIGURE_DEPENDS "*.cu")

# Convert relative paths to absolute paths so the parent can find them
set(ABS_GPU_SOURCES "")
foreach(_file ${GPU_SOURCES})
    get_filename_component(_abs_file "${_file}" ABSOLUTE)
    list(APPEND ABS_GPU_SOURCES "${_abs_file}")
endforeach()

# Send the list up to the numerics/CMakeLists.txt scope
set(CUDA_FILES ${ABS_GPU_SOURCES} PARENT_SCOPE)
