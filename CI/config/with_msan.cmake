# from default, with sanitizer
# no python with clang + sanitizer
include(without_python)
set_option(USE_SANITIZER msan)
set_option(USE_LIBCXX "/libcxx_msan")
set_option(tests_timeout 48)

