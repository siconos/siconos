if(MSVC)
	include(PlatForm/Windows)
	if(CMAKE_SYSTEM_NAME MATCHES "Windows")
		STRING(REGEX REPLACE "\\\\" "/" ENV_PATH "$ENV{PATH}")
		STRING(REGEX REPLACE "\;" "\\\;" ENV_PATH "${ENV_PATH}")
		STRING(REGEX REPLACE "\"" "" ENV_PATH "${ENV_PATH}")
	endif()
        set(EXE_EXT ".exe")

        # Microsoft has its own definition of secure functions ...
        add_compile_options(/D)
        add_compile_options(_CRT_SECURE_NO_WARNINGS)
endif()
