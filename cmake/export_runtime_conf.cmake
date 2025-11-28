#  Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.
#
# Copyright 2025 INRIA.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# --
#[=======================================================================[.rst:

export_runtime_conf
-------------------


 Generate a txt file which contains compile/link info
 to be used at runtime, in the form
 
  INCLUDE_DIRS=...
  DEFINES=...
  FLAGS=...
  LIBS=...
  LINK_DIRS=...

  The resulting file is useful to generate a dedicated header for cppimport when 
  c++ functions are required at runtime, in Python, to run siconos

  Usage
 
  export_runtime_conf(target_name)


#]=======================================================================]
# For a given target, retrieve the real name of the lib and the lib directory.
# Try to handle all possible cases (IMPORTED, local, header-only ...) and apply a specific
# treatment.
# Maybe some are missing ...
#
# Result: set the two following vars in parent scope (caller):
# 
# LIBS_OUT: lib name
# LINK_DIRS_OUT: path where to find the lib
# INC_DIRS_OUT: include dirs
# FLAGS_OUT: compile flags
# DEFS_OUT: compile definitions (-D...)
# 
# 
function(resolve_link_libraries LIBS_OUT LINK_DIRS_OUT INC_DIRS_OUT FLAGS_OUT DEFS_OUT target_name)
    if(TARGET "${target_name}")
        get_target_property(raw_libs ${target_name} INTERFACE_LINK_LIBRARIES)
        get_target_property(raw_incs ${target_name} INTERFACE_INCLUDE_DIRECTORIES)
        get_target_property(raw_opts ${target_name} INTERFACE_COMPILE_OPTIONS)
        get_target_property(raw_defs ${target_name} INTERFACE_COMPILE_DEFINITIONS)
    else()
        return()
    endif()

    set(_libs "") # List of libraries
    set(_dirs "") # List of path to libraries
    set(_inc_dirs "") # List of include dirs
    set(_compile_opts "") # List of compile options
    set(_compile_flags "") # List of compile definitions
 
    if(raw_incs)
        set(_inc_dirs "${raw_incs}")
    endif()
    if(raw_opts)
        set(_compile_flags "${raw_opts}")
    endif()
    if(raw_defs)
        set(_compile_defs "${raw_defs}")
    endif()
    if(NOT raw_libs)
        return()
    endif() 

    foreach(item IN LISTS raw_libs)
        # Skip generator expressions
        if(item MATCHES "^\\$<")
            continue()
        endif()
        # Case 1: item is a known CMake target
        if(TARGET "${item}")
            # Retrieve the target type (STATIC, SHARED, INTERFACE, IMPORTED...)
            get_target_property(type "${item}" TYPE)
            # INTERFACE libraries are header-only 
            if(type STREQUAL "INTERFACE_LIBRARY")
                get_target_property(inter_lib_inc "${item}" INTERFACE_INCLUDE_DIRECTORIES)
                if(inter_lib_inc)
                    list(APPEND _inc_dirs "${inter_lib_inc}")
                endif()
                continue()  
            endif()

            # Retrieve OUTPUT_NAME (if not set, fallback to target name)
            get_target_property(real_name "${item}" OUTPUT_NAME)
            if(NOT real_name)
                set(real_name "${item}")
            endif()
            # Add library name (logical or OUTPUT_NAME)
            # Get actual library file
            get_target_property(imported_loc "${item}" IMPORTED_LOCATION)
            get_target_property(lib_inc "${item}" INTERFACE_INCLUDE_DIRECTORIES)
            if(lib_inc)
                list(APPEND _inc_dirs "${lib_inc}")
            endif() 

            if(imported_loc)
                # For imported targets, IMPORTED_LOCATION gives the full path
                get_filename_component(dir "${imported_loc}" DIRECTORY)
                list(APPEND _dirs "${dir}")
                list(APPEND _libs "${real_name}")
                continue()
                # Rq: Boost::xx are a problem. Imported_location is not defined
                # even if libboost_xx exists
                
            endif()
            if(item MATCHES "::")
                get_target_property(location "${item}" LOCATION)
                if(location)
                    get_filename_component(real_name2 ${location} NAME_WE)
                    get_filename_component(lib_dir ${location} DIRECTORY)
                    string(REGEX REPLACE "^lib" "" real_name2 ${real_name2})
                    list(APPEND _libs "${real_name2}")
                    list(APPEND _dirs "${lib_dir}")
                    continue()
                endif()
                list(APPEND _libs "$<TARGET_FILE_BASE_NAME:${item}>")
            else()
                list(APPEND _libs "${real_name}")
            endif()

            #Local non-imported library
            list(APPEND _dirs "$<TARGET_FILE_DIR:${item}>")
            continue()

        endif()

        # Case: absolute path to a real file
        if(EXISTS "${item}")
            get_filename_component(dir  "${item}" DIRECTORY)
            get_filename_component(name "${item}" NAME_WE)

            list(APPEND _libs "${name}")
            list(APPEND _dirs "${dir}")
            continue()
        endif()

        # Rq: unknown entry -> ignore
    endforeach()
     # Second run to collect properties of deps of input lib
    foreach(subitem IN LISTS raw_libs)
        resolve_link_libraries(sublibs subdirs sub_inc_dirs sub_flags sub_defs "${subitem}")
        list(APPEND _libs "${sublibs}")
        list(APPEND _dirs "${subdirs}")
        list(APPEND _inc_dirs "${sub_inc_dirs}")
        list(APPEND _compile_flags "${sub_flags}")
        list(APPEND _compile_defs "${sub_defs}")
    endforeach()

    # Remove duplicates
    list(REMOVE_DUPLICATES _libs)
    list(REMOVE_DUPLICATES _dirs)
    list(REMOVE_DUPLICATES _inc_dirs)
    list(REMOVE_DUPLICATES _compile_flags)
    list(REMOVE_DUPLICATES _compile_defs)

    # Return results to the caller
    set(${LIBS_OUT}      "${_libs}" PARENT_SCOPE)
    set(${LINK_DIRS_OUT} "${_dirs}" PARENT_SCOPE)
    set(${INC_DIRS_OUT} "${_inc_dirs}" PARENT_SCOPE)
    set(${FLAGS_OUT} "${_compile_flags}" PARENT_SCOPE)
    set(${DEFS_OUT} "${_compile_defs}" PARENT_SCOPE)
endfunction()

function(export_build_conf target_name)

    if(NOT TARGET ${target_name})
        message(FATAL_ERROR "Target ${target_name} does not exist")
    endif()
    resolve_link_libraries(RESOLVED_LIBS RESOLVED_LINK_DIRS RESOLVED_INC_DIRS RESOLVED_FLAGS RESOLVED_DEFS "${target_name}")
    
    set(compile_flags "$<JOIN:$<TARGET_PROPERTY:${target_name},INTERFACE_COMPILE_OPTIONS>,$<SEMICOLON>>")
    list(APPEND compile_flags "${RESOLVED_FLAGS}")
    if(CMAKE_CXX_STANDARD)
        list(APPEND compile_flags "-std=gnu++${CMAKE_CXX_STANDARD}")
    endif()

    set(compile_defs "$<JOIN:$<TARGET_PROPERTY:${target_name},INTERFACE_COMPILE_DEFINITIONS>,$<SEMICOLON>>")
    list(APPEND compile_defs "${RESOLVED_DEFS}")

    set(all_incs "$<JOIN:$<TARGET_PROPERTY:${target_name},INTERFACE_INCLUDE_DIRECTORIES>,$<SEMICOLON>>")
    list(APPEND all_incs "${RESOLVED_INC_DIRS}")

    list(REMOVE_DUPLICATES all_incs)
    list(REMOVE_DUPLICATES compile_defs)
    list(REMOVE_DUPLICATES compile_flags)

    set(out_file "${CMAKE_CURRENT_BINARY_DIR}/cppimport_build_info.txt")
    file(GENERATE
      OUTPUT "${out_file}"
      CONTENT 
      "INCLUDE_DIRS=${all_incs}\n
    # DEFINES=$<JOIN:$<TARGET_PROPERTY:${target_name},INTERFACE_COMPILE_DEFINITIONS>,$<SEMICOLON>>\n
    DEFINES=${compile_defs}\n
    FLAGS=${compile_flags}\n
    LIBS=$<JOIN:${RESOLVED_LIBS},$<SEMICOLON>>\n
    LINK_DIRS=$<JOIN:${RESOLVED_LINK_DIRS},$<SEMICOLON>>\n"
    )
    string(REPLACE "::" "_" target_safe "${target_name}")
    add_custom_target(
        ${target_safe}_generate_cppimport_info
        DEPENDS ${out_file})   
    message(STATUS "Generated runtime config for ${target_name} at ${out_file}")

endfunction()