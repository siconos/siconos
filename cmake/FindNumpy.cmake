# ----------------------------------------------
# Search for numpy headers (arrayobject.h etc .) 
#
# Warning : must be called after FindPythonFull
# ---------------------------------------------
# This module defines
#  NUMPY_INCLUDE_DIR, where to find numpy/arrayobject.h, etc.
#  NUMPY_FOUND, If false, do not try to use numpy headers.
include(LibFindMacros)

if(NUMPY_INCLUDE_DIR)
  # in cache already
  set (NUMPY_FIND_QUIETLY TRUE)
endif (NUMPY_INCLUDE_DIR)

IF(PYTHON_EXECUTABLE)
  EXEC_PROGRAM ("${PYTHON_EXECUTABLE}"
    ARGS -c " \"import numpy; print(numpy.get_include()) \" "
    OUTPUT_VARIABLE NUMPY_INCLUDE_DIR
    RETURN_VALUE NUMPY_NOT_FOUND)
  
  if (NUMPY_INCLUDE_DIR)
    STRING(REGEX REPLACE "\n.*$" "" NUMPY_INCLUDE_DIR "${NUMPY_INCLUDE_DIR}")
    IF(CMAKE_SYSTEM_NAME MATCHES Windows)
      STRING(REGEX REPLACE "\\\\" "/" NUMPY_INCLUDE_DIR "${NUMPY_INCLUDE_DIR}")
    ENDIF()
    set (NUMPY_INCLUDE_DIR ${NUMPY_INCLUDE_DIR} CACHE STRING "Numpy include path")
  endif()
  if (NOT NUMPY_FIND_QUIETLY)
    message (STATUS "Numpy headers found : ${NUMPY_INCLUDE_DIR}")
  endif (NOT NUMPY_FIND_QUIETLY)
endif()

# Final check :
set(Numpy_PROCESS_INCLUDES NUMPY_INCLUDE_DIR)
libfind_process(Numpy)
