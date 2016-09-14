include(CMakeParseArguments)
macro(add_frama_c_test)

  set(options)
  set(oneValueArgs ENTRY)
  set(multiValueArgs INCLUDES)
  cmake_parse_arguments(add_frama_c_test "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  set(SRC ${add_frama_c_test_UNPARSED_ARGUMENTS})

  if(add_frama_c_test_ENTRY)
    set(ENTRY ${add_frama_c_test_ENTRY})
    set(ENTRY_OPT -main=${add_frama_c_test_ENTRY} -lib-entry)
  else()
    set(ENTRY main)
    set(ENTRY_OPT)
  endif()

  set(FRAMA_C_ARGS)
  foreach(_I ${add_frama_c_test_INCLUDES})
    set(FRAMA_C_ARGS ${FRAMA_C_ARGS} -cpp-extra-args=-I${_I})
  endforeach()

  if(NOT FRAMA_C_COMMAND)
    find_program(FRAMA_C_COMMAND frama-c)
  endif()
  
  if(NOT FRAMA_C_PATH)
    execute_process(COMMAND ${FRAMA_C_COMMAND} -print-path OUTPUT_VARIABLE FRAMA_C_PATH)
  endif()

  get_filename_component(SRC_NAME ${SRC} NAME_WE)
  
  set(FRAMA_C_ARGS ${FRAMA_C_ARGS} ${ENTRY_OPT} -cpp-extra-args=-DFUNCODEGEN_CHECK -cpp-extra-args=-I${CMAKE_BINARY_DIR} -cpp-extra-args=-I${CMAKE_CURRENT_SOURCE_DIR} -cpp-extra-args=-I${FRAMA_C_PATH}/libc -cpp-extra-args=-I${FRAMA_C_PATH} -cpp-extra-args=-std=c99 -kernel-msg-key pp -val ${SRC} -val-subdivide-non-linear 1200 -slevel 4096 -print -ocode ${SRC_NAME}_cil.c -then -wp ${SRC} -no-print -wp-model float+int  -wp-timeout 40 -wp-alt-ergo-opt="-backward-compat"  -wp-out temp -wp-verbose 2 -wp-prover alt-ergo  -then -wp -wp-prover z3 -then -wp -wp-prover cvc3 -then -wp -wp-prover cvc4 -then -wp -wp-prover zenon -then -report -then -werror -werror-no-external)

  add_custom_target(
    ${SRC_NAME}-${ENTRY}-frama-c-report
    COMMENT "Static check ${SRC}"
    COMMAND ${FRAMA_C_COMMAND} ${FRAMA_C_ARGS})

  add_test(${SRC_NAME}-${ENTRY}-frama-c-report ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target ${SRC_NAME}-${ENTRY}-frama-c-report)


endmacro()
