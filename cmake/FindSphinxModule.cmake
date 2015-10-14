# Search for a particular sphinx module.
#
# Usage
# -----
# find_sphinx_module(parent module_name REQUIRED)
#
# parent == sphinx or sphinxcontrib
# module_name == name of the required module.
#
# Example:
#
# find_sphinx_module(sphinxcontrib bibtex)
#
# will try in python :
# from sphinxcontrib import bibtex
#
function(find_sphinx_module parent module)
	string(TOUPPER ${module} module_upper)
	if(ARGC GREATER 1 AND ARGV1 STREQUAL "REQUIRED")
	  set(${module}_FIND_REQUIRED TRUE)
	endif()
	execute_process(COMMAND ${PYTHON_EXECUTABLE} -c
	  "import re; from ${parent} import ${module}; print re.compile('/__init__.py.*').sub('',${module}.__file__)"
	  RESULT_VARIABLE _${module}_status
	  OUTPUT_VARIABLE _${module}_location
	  ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)

    	if(NOT _${module}_status)
	  set(python_${module_upper} ${_${module}_location} CACHE STRING
	    "Location of Python module ${module}")
	endif(NOT _${module}_status)

	find_package_handle_standard_args(${module} DEFAULT_MSG _${module}_location)
	set(${module}_FOUND ${${module_upper}_FOUND} PARENT_SCOPE)
endfunction()
