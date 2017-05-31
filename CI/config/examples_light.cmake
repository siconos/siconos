# default configuration for examples only with external, numerics and kernel components
set_components(externals;numerics;kernel)
set_option(DEV_MODE OFF)
set_option(WITH_PYTHON_WRAPPER ON)
set_option(WITH_TESTING OFF)  # We don't want to build/run siconos tests, only examples!


