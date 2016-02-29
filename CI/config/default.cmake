# default configuration for testing

IF(NOT __DEFAULT_CMAKE__)
  SET(__DEFAULT_CMAKE__ TRUE)
  set_option(DEV_MODE ON)
  set_option(WITH_TESTING ON)
  set_option(WITH_SYSTEM_INFO ON)
  set_option(WITH_numerics_TESTING ON)
  set_option(WITH_kernel_TESTING ON)
  set_option(WITH_control_TESTING ON)
  set_option(WITH_mechanics_TESTING ON)
  set_option(WITH_io_TESTING ON)
ENDIF(NOT __DEFAULT_CMAKE__)
