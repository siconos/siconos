#================================================================
# cmake utilities to build external component
#================================================================

set(COMPONENT externals)
message("-- Set up for ${PROJECT_NAME}_${COMPONENT} library ...\n")
# Name used for the target to be built
set(COMPONENT_LIBRARY_NAME ${PROJECT_NAME}_${COMPONENT})

# Define a Fortran interface file (FCMangle.hpp)
#FortranCInterface_VERIFY(QUIET)
if(NOT FortranCInterface_VERIFIED_C)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/DefaultMangle.h.in ${CMAKE_CURRENT_BINARY_DIR}/FCMangle.h)
else()
  FortranCInterface_HEADER(${CMAKE_CURRENT_BINARY_DIR}/FCMangle.h
    MACRO_NAMESPACE "myF2C" 
    SYMBOL_NAMESPACE "myF2C"
    )
endif()

# --- Sources and headers ---
# source directories for numerics
set(${COMPONENT}_DIRS
  .
  blas_lapack
  hairer
  netlib/dftemplates
  netlib/odepack
  optim_misc
  optim_misc/ql0001
  SOL/lumod-c
  sort
  SuiteSparse/CSparse
  SuiteSparse/CSparse_extensions
  tools
)

if(WITH_CXX)
  list(APPEND ${COMPONENT}_DIRS numeric_bindings boost_contribs)
endif()


# Unstable sources. Will be included only
# if WITH_${COMPONENT}_UNSTABLE is true.
set(${COMPONENT}_Unstable_SRCS)

# List of directories of headers not to be installed
set(${COMPONENT}_HDRS_EXCLUDE_DIR PATH_SDK/include SOL/lumod-c)
if(NOT INSTALL_EXTERNALS_HEADERS)
  list(APPEND ${COMPONENT}_HDRS_EXCLUDE_DIR
    blas_lapack
    hairer
    netlib/dftemplates
    netlib/odepack
    optim_misc
    optim_misc/ql0001
    sort
    tools)
endif()
if(NOT INSTALL_EXTERNALS_HEADERS_CSPARSE)
  list(APPEND ${COMPONENT}_HDRS_EXCLUDE_DIR
    SuiteSparse/CSparse
    SuiteSparse/CSparse_extensions)
endif()

# List of directories for which no doxygen doc will be generated
# By default all directories matching "test" are excluded.
set(${COMPONENT}_EXCLUDE_DOXY
  netlib
  hairer
  PATH_SDK/include
  tools)


SET(${COMPONENT}_EXCLUDE_SRCS "SOL/lumod-c/sparselib.c" "SOL/lumod-c/lumod_sparse.c")

set_source_files_properties(SOL/lumod-c/lumod_dense.c PROPERTIES COMPILE_FLAGS -Wno-unused-variable)
set_source_files_properties(SOL/lumod-c/commonlib.c PROPERTIES COMPILE_FLAGS -Wno-unused-variable)
set_source_files_properties(SuiteSparse/CSparse/csparse.c PROPERTIES COMPILE_FLAGS -Wno-unused)



if(NOT HAS_ATLAS_LAPACK)
  set(${COMPONENT}_EXCLUDE_SRCS ${${COMPONENT}_EXCLUDE_SRCS} "blas_lapack/clapack_dtrtrs.c")
endif()

# ---- Final setup for the library ----
# --- set linked libraries and linker language ---
set(${COMPONENT}_LINKER_LANGUAGE C)
# but CXX may be needed with MlcpSimplex
if(MLCPSIMPLEX_FOUND) 
  set(${COMPONENT}_LINKER_LANGUAGE CXX)
endif()
set(${COMPONENT}_LINK_LIBRARIES ${SICONOS_LINK_LIBRARIES})

# Windows stuff
include(WindowsExternalsSetup)

if(BUILD_AS_CPP)
  file(GLOB_RECURSE C_FILES ${CMAKE_CURRENT_SOURCE_DIR} *.c)
  set_source_files_properties(${C_FILES} PROPERTIES LANGUAGE CXX)
  set(${COMPONENT}_LINKER_LANGUAGE CXX)
endif()

include(LibraryProjectSetup)
library_project_setup()



# --- specific install ----
#install(FILES ${CMAKE_CURRENT_BINARY_DIR}/FCMangle.h DESTINATION include/${PROJECT_NAME})

if(HAVE_SICONOS_MECHANICS AND WITH_RENDERER)
  # This files is not installed as a program, even though it is one.
  # Remove the comment if you know how it should be installed --xhub
  #  configure_file(io/SimpleGui.py ${SICONOS_SWIG_ROOT_DIR}/io/SimpleGui.py @ONLY) 


  file(GLOB rendererFiles RELATIVE ${CMAKE_SOURCE_DIR}/externals/renderer  ${CMAKE_SOURCE_DIR}/externals/renderer/img/*.*)
  foreach(rendererFile ${rendererFiles})
    set(srcRendererPath ${CMAKE_SOURCE_DIR}/externals/renderer/${rendererFile})
    if(NOT IS_DIRECTORY ${srcRendererPath})     
      install(FILES  ${srcRendererPath} DESTINATION share/siconos/renderer/img)
    endif(NOT IS_DIRECTORY ${srcRendererPath})
  endforeach(rendererFile)

  file(GLOB rendererFiles RELATIVE ${CMAKE_SOURCE_DIR}/externals/renderer/  ${CMAKE_SOURCE_DIR}/externals/renderer/threeJS_libraries/*.*)
  foreach(rendererFile ${rendererFiles})
    set(srcRendererPath ${CMAKE_SOURCE_DIR}/externals/renderer/${rendererFile})
    if(NOT IS_DIRECTORY ${srcRendererPath})     
      install(FILES  ${srcRendererPath} DESTINATION share/siconos/renderer/threeJS_libraries  )
    endif(NOT IS_DIRECTORY ${srcRendererPath})
  endforeach(rendererFile)

  if(INSTALL_PYTHON_SYMLINKS)
    message("Setting up symlink install targets for externals Python executables")
    install(CODE "execute_process(COMMAND sh -c \"test -e '${CMAKE_INSTALL_PREFIX}/bin/siconos_renderer' || ln -vs '${CMAKE_CURRENT_SOURCE_DIR}/renderer/renderer.py' '${CMAKE_INSTALL_PREFIX}/bin/siconos_renderer'\")")
  else()
    install(PROGRAMS renderer/renderer.py
      DESTINATION bin RENAME siconos_renderer)
  endif()

endif()

if(WITH_SERIALIZATION)
  install(FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/boost_contribs/boost/serialization/unordered_collections_load_imp.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/boost_contribs/boost/serialization/unordered_collections_save_imp.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/boost_contribs/boost/serialization/unordered_map.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/boost_contribs/boost/serialization/unordered_set.hpp
    DESTINATION include/${PROJECT_NAME}/boost/serialization)
endif()




# --- tests ---
include(${COMPONENT}_tests.cmake)
