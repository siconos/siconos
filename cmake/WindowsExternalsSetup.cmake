IF(MSVC AND "${CMAKE_Fortran_COMPILER}" MATCHES "gfortran")
  SET(GCC_A "")
  EXECUTE_PROCESS(COMMAND ${CMAKE_Fortran_COMPILER} -print-file-name=libgfortran.dll.a
	  OUTPUT_VARIABLE LIBGFORTRAN_DLL_A)
  find_program(CYGPATH cygpath)
  IF(CYGPATH)
    EXECUTE_PROCESS(COMMAND cygpath -m ${LIBGFORTRAN_DLL_A}
      OUTPUT_VARIABLE LIBGFORTRAN)
  ELSE(CYGPATH) # Hope for the best
    SET(LIBGFORTRAN "${LIBGFORTRAN_DLL_A}")
  ENDIF(CYGPATH)
  STRING(REGEX REPLACE "\n.*$" "" LIBGFORTRAN ${LIBGFORTRAN})
  MESSAGE(STATUS "libgfortran.dll.a :: ${LIBGFORTRAN}")
  SET(${COMPONENT}_LINK_LIBRARIES ${${COMPONENT}_LINK_LIBRARIES} ${LIBGFORTRAN})
  GET_FILENAME_COMPONENT(GFORTRAN_DIR ${CMAKE_Fortran_COMPILER} PATH)
  IF(NOT CMAKE_AR)
    SET(CMAKE_AR "${GFORTRAN_DIR}/ar")
  ENDIF(NOT CMAKE_AR)
  target_compile_options(externals PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:static>)
  target_compile_options(externals PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-static-libgcc>)
  target_compile_options(externals PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-static-libgfortran>)
  target_compile_options(externals PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-fno-stack-check>)
  target_compile_options(externals PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-fno-stack-protector>)
  target_compile_options(externals PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-mno-stack-arg-probe>)
   # XXX No test :( -- xhub
endif()

