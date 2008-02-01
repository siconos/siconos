ENABLE_TESTING()

ADD_EXECUTABLE(MLCPtest main_mlcp.c)

# mic mac
TARGET_LINK_LIBRARIES(MLCPtest SiconosNumerics ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} gfortran lapack_atlas)

ADD_TEST(T1 MLCPtest)

