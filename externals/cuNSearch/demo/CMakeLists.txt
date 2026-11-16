add_executable(Demo
    main.cu
)
add_dependencies(Demo cuNSearch)
target_link_libraries(Demo PRIVATE cuNSearch)
