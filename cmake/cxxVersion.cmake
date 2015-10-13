# Based on the Qt 5 processor detection code, so should be very accurate
# https://qt.gitorious.org/qt/qtbase/blobs/master/src/corelib/global/qprocessordetection.h
# Currently handles arm (v5, v6, v7), x86 (32/64), ia64, and ppc (32/64)

# Regarding POWER/PowerPC, just as is noted in the Qt source,
# "There are many more known variants/revisions that we do not handle/detect."

set(cxx_detect_code "
#if __cplusplus >= 201103L
#error CXXVERSION 201103L
#elif __cplusplus >= 199711L
#error CXXVERSION 199711L
#else
#error CXXVERSION 000000L
#endif
")

# Set ppc_support to TRUE before including this file or ppc and ppc64
# will be treated as invalid architectures since they are no longer supported by Apple

function(detect_cxx_version output_var)
        file(WRITE "${CMAKE_BINARY_DIR}/cxxversion.cpp" "${cxx_detect_code}")

        enable_language(CXX)

        try_run(
            run_result_unused
            compile_result_unused
            "${CMAKE_BINARY_DIR}"
            "${CMAKE_BINARY_DIR}/cxxversion.cpp"
            COMPILE_OUTPUT_VARIABLE CXXVERSION
        )

        # Parse the architecture name from the compiler output
        string(REGEX MATCH "CXXVERSION ([a-zA-Z0-9_]+)" CXXVERSION "${CXXVERSION}")

        # Get rid of the value marker leaving just the architecture name
        string(REPLACE "CXXVERSION " "" CXXVERSION "${CXXVERSION}")

        set(${output_var} "${CXXVERSION}" CACHE STRING "C++ standart used to compile the process" FORCE)
endfunction()
