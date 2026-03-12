

add_definitions(-D_USE_MATH_DEFINES)

include_directories("../include")

add_executable(Demo
    main.cpp
)

target_link_libraries(Demo CompactNSearch)
