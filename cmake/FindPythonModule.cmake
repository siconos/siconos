# Search for a Python module.
#
# Usage :
#
#  find_python_module(mpi4py REQUIRED)
#  find_python_module(sphinx)
#  find_python_module(mpi4py INCLUDE)
# 
#  Warning : use ${PYTHON_EXECUTABLE} as python interpreter
#
#  If INCLUDE options is provided, it means that the function
#  is supposed to check for the existence of <path-to-module>/include
#  and set ${module}_INCLUDE_DIR cache variable.


include(FindPackageHandleStandardArgs)
function(find_python_module module)
  set(options REQUIRED INCLUDES) # If INCLUDE options is provided
  set(oneValueArgs VERSION)
  cmake_parse_arguments(${module}_FIND "${options}" "${oneValueArgs}" "" ${ARGN} )
  execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "import ${module} as name; print(name.__file__)"
    RESULT_VARIABLE ${module}_FIND_RESULT     # Return code from command above
    OUTPUT_VARIABLE ${module}_FIND_OUTPUT     # Standard output form command above
    ERROR_QUIET # Ignores quietly standard error
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )
  
  if(NOT ${module}_FIND_RESULT) # Return code == 0 means that things have gone well
    set(${module}_file ${${module}_FIND_OUTPUT} CACHE STRING "Python ${module} module file.")
  else()
    set(${module}_file ${module}-NOTFOUND)
  endif()

  # Save version
  execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "import ${module} as name; print(name.__version__)"
    RESULT_VARIABLE ${module}_FIND_RESULT     # Return code from command above
    OUTPUT_VARIABLE ${module}_VERSION    # Standard output form command above
    ERROR_QUIET # Ignores quietly standard error
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )
  if(${module}_FIND_RESULT)
    execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "import ${module} as name; print(name.VERSION)"
      RESULT_VARIABLE ${module}_FIND_RESULT     # Return code from command above
      OUTPUT_VARIABLE ${module}_VERSION    # Standard output form command above
      ERROR_QUIET # Ignores quietly standard error
      OUTPUT_STRIP_TRAILING_WHITESPACE
      )
  endif()
  
  set(${module}_REQUIRED_VARS ${module}_file ${module}_PYTHONPATH)
    
  if(${module}_file)
    # module directory
    get_filename_component(${module}_path ${${module}_file} DIRECTORY)
    # path to the module (i.e. what must be added to PYTHONPATH or sys.path)
    get_filename_component(${module}_PYTHONPATH ${${module}_path} DIRECTORY)
    
    if(${module}_FIND_INCLUDES)
      find_file(${module}_INCLUDE_DIR include PATHS ${${module}_path} NO_DEFAULT_PATH)
      list(APPEND ${module}_REQUIRED_VARS ${module}_INCLUDE_DIR)
    endif()
  endif()
  find_package_handle_standard_args(${module}
    REQUIRED_VARS ${${module}_REQUIRED_VARS}
    VERSION_VAR ${module}_VERSION)
 
  if(${module}_FOUND)
    message("-- Found python package ${${module}_file}, version ${${module}_VERSION}")
  endif()
endfunction(find_python_module)
