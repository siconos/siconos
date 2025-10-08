# Generate a txt file which contains compile/link info
# to be used at runtime, in the form
# 
#  INCLUDE_DIRS=...
#  DEFINES=...
#  FLAGS=...
#  LIBS=...
#  LINK_DIRS=...
#
#  The resulting file is useful to generate a dedicated header for cppimport when 
#  c++ functions are required at runtime, in Python, to run siconos
#
#  Usage
# 
#  export_runtime_conf()
#
function(export_runtime_conf target)
    set(oneValueArgs VAR) # output variable name
    set(multiValueArgs DIRS EXTS)
    set(options RECURSIVE)
    cmake_parse_arguments(collect "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
    
    # Look for target properties
    #get_target_property(_inc_dirs ${target} INCLUDE_DIRECTORIES)
    #get_target_property(_compile_defs ${target} COMPILE_DEFINITIONS)
    #get_target_property(_compile_opts ${target} COMPILE_OPTIONS)
    get_target_property(_link_libs  ${target} LINK_LIBRARIES)
    #get_target_property(_link_dirs  ${target} LINK_DIRECTORIES)



    if (CMAKE_CXX_STANDARD)
        list(APPEND _compile_opts "-std=gnu++${CMAKE_CXX_STANDARD}")
    endif()

    # clean lists
    foreach(var _inc_dirs _compile_defs _compile_opts _link_libs _link_dirs)
        list(FILTER ${var} EXCLUDE REGEX "NOTFOUND")
        list(REMOVE_DUPLICATES ${var})
    endforeach()

    # Export in a file that can be processed by python
    set(_outfile "${CMAKE_BINARY_DIR}/siconos_runtime_build_info.txt")
    file(WRITE ${_outfile} "INCLUDE_DIRS=${_inc_dirs}\n")
    file(APPEND ${_outfile} "DEFINES=${_compile_defs}\n")
    file(APPEND ${_outfile} "FLAGS=${_compile_opts}\n")
    file(APPEND ${_outfile} "LIBS=${_link_libs}\n")
    file(APPEND ${_outfile} "LINK_DIRS=${_link_dirs}\n")

    message(STATUS "✅ Siconos runtime build info has been written in: ${_outfile}")

endfunction()


# For a given target, collect info regarding 
# - compile options and definitions
# - linked libraries
# - include path ...
# 
# Used to prepare a txt file for later cppimport 
function(export_properties)

    set(oneValueArgs VAR) 
    cmake_parse_arguments(lib "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

    set(_inc_dirs)
    set(lib_defs)
    set(lib_opts)
    set(_link_libs)
    set(_resolved_libs)
    set(lib_location)
    
    # includes - Set _inc_dirs
    get_target_property(lib_inc ${lib} INTERFACE_INCLUDE_DIRECTORIES)
    if(lib_inc)
      list(APPEND _inc_dirs ${lib_inc})
    else()
      get_target_property(lib_inc ${lib} IMPORTED_LOCATION)
      if(lib_inc)
        get_filename_component(lib_dir ${lib_inc} DIRECTORY)
        list(APPEND _inc_dirs ${lib_dir})
        message("Imported special case inc : ${lib_dir}")
      endif()
    endif()

    # compile options/def
    get_target_property(lib_defs ${lib} INTERFACE_COMPILE_DEFINITIONS)
    get_target_property(lib_opts ${lib} INTERFACE_COMPILE_OPTIONS)
  
    # library location, name ...
    get_target_property(lib_dirs ${lib} INTERFACE_LINK_DIRECTORIES)
    list(APPEND _link_dirs ${lib_dirs})
    get_target_property(lib_link ${lib} INTERFACE_LINK_LIBRARIES)
    get_target_property(lib_imp_inc ${lib} IMPORTED_IMPLIB)
    list(APPEND _link_libs ${lib_link})
    list(APPEND _link_libs ${lib_imp_inc})


    get_target_property(lib_imported_inc ${lib} IMPORTED_LOCATION)
    get_target_property(lib_location ${lib} LOCATION)

    if(lib_location)
        get_filename_component(real_name ${lib_location} NAME_WE)
        get_filename_component(lib_dir ${lib_location} DIRECTORY)
        string(REGEX REPLACE "^lib" "" real_name ${real_name})  # enleve le prefixe lib
        list(APPEND _resolved_libs ${real_name})
        list(APPEND _link_dirs ${lib_dir})
    else()
        list(APPEND _resolved_libs ${lib})      
    endif()
        
    set(${lib}_inc_dirs ${_inc_dirs} PARENT_SCOPE)
    set(${lib}_compile_defs ${lib_defs} PARENT_SCOPE)
    set(${lib}_compile_opts ${lib_opts} PARENT_SCOPE)
    set(${lib}_link_dirs ${_link_dirs} PARENT_SCOPE)
    set(${lib}_link_libs ${_link_libs} PARENT_SCOPE)
    set(${lib}_resolved_libs ${_resolved_libs} PARENT_SCOPE)
    set(${lib}_real_name ${real_name} PARENT_SCOPE)
    set(${lib}_lib_location ${lib_location} PARENT_SCOPE)  
endfunction()

function (export_runtime_conf MyTarget)

    # First run to collect input lib properties
    get_target_property(_link_libs  ${MyTarget} LINK_LIBRARIES)

    foreach(lib ${_link_libs})
        export_properties(VAR ${lib})
        list(APPEND all_inc_dirs ${${lib}_inc_dirs} )
        list(APPEND all_compile_defs ${${lib}_compile_defs} )
        list(APPEND all_compile_opts ${${lib}_compile_opts} )
        list(APPEND all_link_dirs ${${lib}_link_dirs} )
        list(APPEND all_link_libs ${${lib}_link_libs} )
        list(APPEND all_resolved_libs ${${lib}_resolved_libs} )
        list(APPEND all_real_name ${${lib}_real_name} )
        list(APPEND all_lib_location ${${lib}_lib_location} )

        foreach(var all_inc_dirs all_compile_opts all_compile_defs all_link_libs
                all_link_dirs all_resolved_libs  all_real_name all_lib_location)
            list(FILTER ${var} EXCLUDE REGEX "NOTFOUND")
            list(REMOVE_DUPLICATES ${var})
        endforeach()
    endforeach()

    # Second run to collect properties of deps of input lib
    foreach(lib IN LISTS all_link_libs)
        list(FIND _link_libs ${lib} index)
        if(index EQUAL -1)
            export_properties(VAR ${lib})
            list(APPEND all_inc_dirs ${${lib}_inc_dirs} )
            list(APPEND all_compile_defs ${${lib}_compile_defs} )
            list(APPEND all_compile_opts ${${lib}_compile_opts} )
            list(APPEND all_link_dirs ${${lib}_link_dirs} )
            list(APPEND all_link_libs ${${lib}_link_libs} )
            list(APPEND all_resolved_libs ${${lib}_resolved_libs} )
            list(APPEND all_real_name ${${lib}_real_name} )
            list(APPEND all_lib_location ${${lib}_location} )
            foreach(var all_inc_dirs all_compile_opts all_compile_defs all_link_libs
                    all_link_dirs all_resolved_libs  all_real_name all_lib_location)
                list(FILTER ${var} EXCLUDE REGEX "NOTFOUND")
                list(REMOVE_DUPLICATES ${var})
            endforeach()
        endif()
    endforeach()
    if (CMAKE_CXX_STANDARD)
        list(APPEND all_compile_opts "-std=gnu++${CMAKE_CXX_STANDARD}")
    endif()

    set(all_link_libs ${all_resolved_libs})

    # Export into a file 
    set(_outfile "${CMAKE_BINARY_DIR}/cppimport_build_info.txt")
    file(WRITE ${_outfile} "INCLUDE_DIRS=${all_inc_dirs}\n")
    file(APPEND ${_outfile} "DEFINES=${all_compile_defs}\n")
    file(APPEND ${_outfile} "FLAGS=${all_compile_opts}\n")
    file(APPEND ${_outfile} "LIBS=${all_link_libs}\n")
    file(APPEND ${_outfile} "LINK_DIRS=${all_link_dirs}\n")

    message(STATUS "✅ Done export to: ${_outfile}")
endfunction()