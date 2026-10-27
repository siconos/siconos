# cmake -DA=<exe> -DB=<exe> -DWD=<dir> -P compare_outputs.cmake
#
# runs A then B in WD: both must succeed and print the same output. The
# relation and index sets oracle: A and B are the same test built without
# and with the sets (VKERNEL_RELATION_SETS, VKERNEL_INDEX_SETS).

cmake_minimum_required(VERSION 3.17)  # foreach(IN ZIP_LISTS), list policies

file(MAKE_DIRECTORY ${WD})

foreach(exe A B)
  execute_process(COMMAND ${${exe}}
    WORKING_DIRECTORY ${WD}
    OUTPUT_VARIABLE out_${exe}
    ERROR_VARIABLE err_${exe}
    RESULT_VARIABLE res_${exe})
  if(NOT res_${exe} EQUAL 0)
    message(FATAL_ERROR "${${exe}} failed (${res_${exe}}):\n${out_${exe}}\n${err_${exe}}")
  endif()
endforeach()

if(NOT out_A STREQUAL out_B)
  string(REPLACE "\n" ";" lines_A "${out_A}")
  string(REPLACE "\n" ";" lines_B "${out_B}")
  list(LENGTH lines_A n_A)
  list(LENGTH lines_B n_B)
  set(line 0)
  foreach(a b IN ZIP_LISTS lines_A lines_B)
    if(NOT a STREQUAL b)
      message(FATAL_ERROR
        "outputs differ at line ${line} (${n_A} / ${n_B} lines):\n"
        "${A}: ${a}\n${B}: ${b}")
    endif()
    math(EXPR line "${line} + 1")
  endforeach()
  message(FATAL_ERROR "outputs differ in length: ${n_A} / ${n_B} lines")
endif()
