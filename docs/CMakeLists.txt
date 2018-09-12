# ============================================================
# Prepare  config to generate siconos documentation :
#
# - automatic doxygen documentation, from source code
# - sphinx documentation for manuals
# 
# ============================================================
include(FindPythonModule)
include(FindSphinxModule)

# Devel option : activate this to generate rst from doxygen xml outputs
# using doxyrest. 
option(USE_DOXYREST "Build rst files from doxygen outputs (tool = doxyrest). Default = OFF" OFF)
# Devel option : activate this to generate rst from doxygen xml outputs
# using exhale
option(USE_EXHALE "Build rst files from doxygen outputs (tool = exhale). Default = OFF" OFF)

# Devel option : generate pure doxygen doc (no sphinx), extracting all and with graphs.
option(USE_DEVEL_DOXYGEN "Build doxygen (html, latex ...) doc, extracting all. Devel purpose. Default = OFF" OFF)

set(DOC_ROOT_DIR ${CMAKE_BINARY_DIR}/docs/build CACHE INTERNAL "Root path to generated documentation.")
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/build/html)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doctools.py ${CMAKE_BINARY_DIR}/share/doctools.py COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/replace_py_escape.sh ${CMAKE_BINARY_DIR}/share/replace_py_escape.sh COPYONLY)

#  ====== doxygen documentation setup ======
# This setup is required if one at least of
# the following is needed:
# - doxygen warnings (see WITH_DOXYGEN_WARNINGS)
# - python docstrings from doxygen (xml) (see WITH_DOXY2SWIG)
# - doxygen documentation generated from source files (see WITH_DOCUMENTATION)
if(USE_DOXYGEN)

  # Search doxygen and set config
  find_package(Doxygen REQUIRED)
  # Output path for doxygen documentation
  if(NOT DOXYGEN_OUTPUT)
    set(DOXYGEN_OUTPUT ${DOC_ROOT_DIR}/doxygen CACHE INTERNAL "Path to doxygen outputs.")
  endif()

  file(MAKE_DIRECTORY ${DOXYGEN_OUTPUT}/man)
  
  # --- Config files and their options. --
  # --> a lot for doxygen warnings (from doxy_warnings.config.in)
  # --> one for documentation (html and xml to sphinx) (from doxy.config.in)

  # == set configuration file for doxygen warnings ==
  #  - Results in binary_dir/doxygen_warnings
  #  - input config from config/doxy_warnings.config
  #  - only xml output.
  if(WITH_DOXYGEN_WARNINGS)
    file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/doxygen_warnings)
    file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/doxygen_warnings/man)
    file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/doxygen_warnings/xml)
    set(DOXY_WARNINGS_CONFIG ${CMAKE_CURRENT_SOURCE_DIR}/config/doxy_warnings.config.in CACHE INTERNAL "config file for doxygen warnings output.")
    # All the configure process will be done for each component
    # thanks to 'include(doxygen_warnings)'. See LibraryProjectSetup.
  endif()

  # == set configuration file for doxygen doc from sources ==
  #  - Results in binary_dir/docs/
  #  - input config from config/doxy.config
  #  - only html output
  # 'inputs' are updated by each component, during call to update_doxy_config_file macro.
  # Doc will be built when 'make doxygen' is called.
  if(WITH_DOCUMENTATION)
    # config file name
    set(DOXY_CONFIG "${CMAKE_CURRENT_BINARY_DIR}/config/doxy.config" CACHE INTERNAL "Doxygen configuration file : used to produce html (doxygen only) and xml files for sphinx (breathe/exhale).")
    # The generation of the config and the creation of the target will be done later,
    # after the update of inputs by each component, with macro 'finalize_doc'
    
  endif()
endif()

  
#  ====== Create a Target to generate sphinx documentation ======
if(USE_SPHINX)
  if(USE_DOXYREST)
    find_package(doxyrest CONFIGS doxyrest_config.cmake REQUIRED)
    if(doxyrest_FOUND)
      set(HAS_DOXYREST "True" CACHE INTERNAL "True if doxyrest required and found.")
      message(STATUS "Doxyrest found : ${DOXYREST_EXE}")
    endif()
  else()
    set(HAS_DOXYREST "False" CACHE INTERNAL "True if doxyrest required and found.")
  endif()
  find_python_module(breathe REQUIRED)
  if(USE_EXHALE)
    find_python_module(exhale REQUIRED)
    set(HAS_EXHALE "True" CACHE INTERNAL "True if exhale required and found.")
  else()
    set(HAS_EXHALE "False" CACHE INTERNAL "True if exhale required and found.")
  endif()
  find_python_module(sphinx REQUIRED)
  find_sphinx_module(sphinxcontrib bibtex REQUIRED)
  #find_python_module(numpydoc REQUIRED)
  #find_sphinx_module(sphinxcontrib doxylink REQUIRED)
  find_sphinx_module(sphinxcontrib youtube)
  if(NOT ${youtube_FOUND})
    # This python package can not be found using standard ways such as conda or pip,
    # so we force download and install.
    message("Can not found sphinxcontrib-youtube. It will be downloaded and installed ...")
    execute_process(
      COMMAND ${PYTHON_EXECUTABLE} -m pip install git+https://github.com/sphinx-contrib/youtube.git)
    find_sphinx_module(sphinxcontrib youtube REQUIRED)
  endif()
  # Sphinx theme for documentation
  #find_python_module(sphinx-bootstrap-theme REQUIRED)
  find_python_module(sphinx_bootstrap_theme REQUIRED)
  configure_file (
    "${CMAKE_CURRENT_SOURCE_DIR}/sphinx/index.rst.in"
    "${CMAKE_CURRENT_BINARY_DIR}/sphinx/index.rst" @ONLY)

  configure_file (
    "${CMAKE_CURRENT_SOURCE_DIR}/sphinx/siconos.bib"
    "${CMAKE_CURRENT_BINARY_DIR}/sphinx/siconos.bib" COPYONLY)

  configure_file(${CMAKE_SOURCE_DIR}/externals/swig/parse_doxygen.py ${CMAKE_BINARY_DIR}/share/parse_doxygen.py COPYONLY)
  
  # List and copy all rst files from source dir to binary dir
  file(GLOB_RECURSE files_rst RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}/sphinx
    ${CMAKE_CURRENT_SOURCE_DIR}/sphinx/*.rst)
  foreach(filepath ${files_rst})
    #get_filename_component(filename ${filepath} NAME)
    
    configure_file (
      "${CMAKE_CURRENT_SOURCE_DIR}/sphinx/${filepath}"
      "${CMAKE_CURRENT_BINARY_DIR}/sphinx/${filepath}" COPYONLY)
  endforeach()
  
  # html menu
  #configure_file (
  #  "${CMAKE_CURRENT_SOURCE_DIR}/sphinx/_templates/menu.html"
  #  "${CMAKE_CURRENT_BINARY_DIR}/sphinx/_templates/menu.html" COPYONLY)

  # css
  configure_file (
    "${CMAKE_CURRENT_SOURCE_DIR}/sphinx/_static/biblio.css"
    "${CMAKE_CURRENT_BINARY_DIR}/sphinx/_static/biblio.css" COPYONLY)

  
  # --- Build a target for sphinx ---
  # --> use make html or make latex ... to build documentation.
  set(SPHINX_EXECUTABLE sphinx-build)
  set(SPHINX_LOG_FILE ${CMAKE_BINARY_DIR}/docs/sphinx_warnings.log)    
  if(HAS_EXHALE)
    set(SPHINX_PARAMETERS
      -b html # html output
      -d build/doctrees # path to doctree files
      ${CMAKE_CURRENT_BINARY_DIR}/sphinx  # path to rst source
      -t WithExhale
      )
    
  else()
    set(SPHINX_PARAMETERS
      -w ${SPHINX_LOG_FILE}
      -b html # html output
      -d build/doctrees # path to doctree files
      ${CMAKE_CURRENT_BINARY_DIR}/sphinx  # path to rst source
      -t WithoutExhale
      )
  endif()

  # -- Target to copy figures required by sphinx from source to binary dir --
  #  make copy-sphinx-figures
  add_custom_target(copy-sphinx-figures
    COMMAND cmake -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/sphinx/figures ${CMAKE_CURRENT_BINARY_DIR}/sphinx/figures
    )
  
  add_custom_target(create_rst
    COMMAND ${CMAKE_COMMAND} -E env PYTHONPATH=${CMAKE_BINARY_DIR}/share ${PYTHON_EXECUTABLE} -c
    "import doctools; doctools.build_api_main('${CMAKE_CURRENT_BINARY_DIR}/sphinx')"
    VERBATIM
    COMMENT "Create main rst files for C/C++ and python API.")
  
  # -- To build html doc : make html --
  add_custom_target(html
    COMMAND ${SPHINX_EXECUTABLE} ${SPHINX_PARAMETERS} 2> errors.txt ${CMAKE_CURRENT_BINARY_DIR}/build/html
    VERBATIM
    DEPENDS copy-sphinx-figures create_rst
    )

  add_custom_command(TARGET html POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E env PYTHONPATH=${CMAKE_BINARY_DIR}/share ${PYTHON_EXECUTABLE} -c
      "import doctools; doctools.filter_sphinx_warnings('${SPHINX_LOG_FILE}')"
      VERBATIM
      COMMENT "Filtering sphinx logfile into ${SPHINX_LOG_FILE}.")

  add_custom_command(TARGET html POST_BUILD
      COMMENT "Done. Sphinx doc has been generated in ${DOC_ROOT_DIR}/html")

  # Option to be tested later : generate docs for a list of tags, branches ...
  add_custom_target(html-multiversion
    COMMAND sphinx-versioning -g ${CMAKE_SOURCE_DIR} build docs/sphinx 
    ${CMAKE_CURRENT_BINARY_DIR}/build/html
    )
  
  # -- To use doxyrest to generate rst files from doxygen (xml) outputs : make doxyrest --
  # Note FP : test purpose only, exhale/breathe seems to be better than doxyrest.
  if(HAS_DOXYREST)
    add_custom_target(doxyrest
      COMMAND ${DOXYREST_EXE} ${DOXYGEN_OUTPUT}/xml/index.xml -o ${CMAKE_CURRENT_BINARY_DIR}/sphinx/from_doxygen/index.rst -F ${DOXYREST_FRAME_DIR} -f c_index.rst.in
      )
    add_dependencies(html doxyrest)
  endif()

  # -- To build latex doc : make latex -- 
  set(SPHINX_LATEX_PARAMETERS
    -b latex # latex output
    -d build/doctrees # path to doctree files
    ${CMAKE_CURRENT_BINARY_DIR}/sphinx  # path to rst source
    )
  
  add_custom_target(latex
    COMMAND ${SPHINX_EXECUTABLE} ${SPHINX_LATEX_PARAMETERS}
    ${CMAKE_CURRENT_BINARY_DIR}/build/latex
    )

  add_custom_target(doc DEPENDS html)
  
  # # -- Target to build python api documention with apidoc --
  # # make apidoc
  # if(WITH_PYTHON_WRAPPER)
  #   set(SPHINX_APIDOC sphinx-apidoc)
  #   set(GENERATED_DOC ${CMAKE_CURRENT_BINARY_DIR}/sphinx/reference)
  #   set(APIDOC_PARAMETERS -e -M -o ${GENERATED_DOC})
  #   # -e : Put documentation for each module on its own page
  #   # -P : Include “_private” modules. (should we??)
  #   # -M : Put module documentation before submodule documentation.
  #   # -a : append documented module path to sys.path
  #   list(GET COMPONENTS -1 _last_comp)
  #   # apidoc must be run AFTER swig process --> depends on components
  #   add_custom_target(apidoc
  #     COMMAND ${CMAKE_COMMAND} -E env "SPHINX_APIDOC_OPTIONS='members','undoc-members'" ${SPHINX_APIDOC} ${APIDOC_PARAMETERS}
  #     ${CMAKE_BINARY_DIR}/wrap/${SICONOS_PYTHON_PACKAGE}
  #     DEPENDS ${_last_comp}
  #     )

  #   # Sphinx html gen. depends on apidoc outputs.
  #   if(WITH_DOXY2SWIG)
  #     #add_dependencies(html apidoc)
  #   endif()
  # endif()

endif()
