# Search for a particular python module.
# Sources : see http://www.cmake.org/pipermail/cmake/2011-January/041666.html
#
# Usage : find_python_module(mpi4py REQUIRED)
#
include(FindPackageHandleStandardArgs)
function(find_python_module module)
	string(TOUPPER ${module} module_upper)
	if(ARGC GREATER 1 AND ARGV1 STREQUAL "REQUIRED")
	  set(${module}_FIND_REQUIRED TRUE)
	endif()
	# A module's location is usually a directory, but for binary modules
	# it's a .so file.
	execute_process(COMMAND ${PYTHON_EXECUTABLE} -c
		"import re, ${module}; print(re.compile('/__init__.py.*').sub('',${module}.__file__))"
		RESULT_VARIABLE _${module}_status
		OUTPUT_VARIABLE _${module}_location
		ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)

	if(NOT _${module}_status)
		set(python_${module_upper} ${_${module}_location} CACHE STRING
			"Location of Python module ${module}")
	endif(NOT _${module}_status)

	execute_process(COMMAND ${PYTHON_EXECUTABLE} -c
		"import ${module}; print(${module}.__version__)"
	  OUTPUT_VARIABLE _${module}_version
	  ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)

	find_package_handle_standard_args(${module} REQUIRED_VARS _${module}_location VERSION_VAR _${module}_version)
	set(${module}_FOUND ${${module_upper}_FOUND} PARENT_SCOPE)
	set(${module}_VERSION ${_${module}_version} PARENT_SCOPE)
endfunction(find_python_module)
