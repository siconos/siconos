# Based on the Qt 5 processor detection code, so should be very accurate
# https://qt.gitorious.org/qt/qtbase/blobs/master/src/corelib/global/qprocessordetection.h
# Currently handles arm (v5, v6, v7), x86 (32/64), ia64, and ppc (32/64)

# Regarding POWER/PowerPC, just as is noted in the Qt source,
# "There are many more known variants/revisions that we do not handle/detect."

set(c_detect_code "
#if __STDC_VERSION__ >= 201112L
#error STDC_VERSION 2011012L
#else
#error STDC_VERSION 000000L
#endif                                                                                                                                                                                                                                      ")

# Set ppc_support to TRUE before including this file or ppc and ppc64
# will be treated as invalid architectures since they are no longer supported by Apple

function(detect_c_version output_var)
        file(WRITE "${CMAKE_BINARY_DIR}/cversion.c" "${c_detect_code}")

        enable_language(C)

        try_run(
            run_result_unused
            compile_result_unused
            "${CMAKE_BINARY_DIR}"
            "${CMAKE_BINARY_DIR}/cversion.c"
            COMPILE_OUTPUT_VARIABLE CVERSION
        )

        # Parse the architecture name from the compiler output
        string(REGEX MATCH "STDC_VERSION ([a-zA-Z0-9_]+)" CVERSION "${CVERSION}")

        # Get rid of the value marker leaving just the architecture name
        string(REPLACE "STDC_VERSION " "" CVERSION "${CVERSION}")

        set(${output_var} "${CVERSION}" PARENT_SCOPE)
endfunction()
