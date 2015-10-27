IF(MSVC)
	include(PlatForm/Windows)
	IF(CMAKE_SYSTEM_NAME MATCHES "Windows")
		STRING(REGEX REPLACE "\\\\" "/" ENV_PATH "$ENV{PATH}")
		STRING(REGEX REPLACE "\;" "\\\;" ENV_PATH "${ENV_PATH}")
		STRING(REGEX REPLACE "\"" "" ENV_PATH "${ENV_PATH}")
	ENDIF()
        SET(EXE_EXT ".exe")

        # Microsoft has its own definition of secure functions ...
        APPEND_CXX_FLAGS("/D _CRT_SECURE_NO_WARNINGS")
ENDIF(MSVC)
