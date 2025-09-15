include(tools4tests)

if(WITH_TESTING)

  add_custom_target(externals-tests echo "Start externals tests")
  
  if(HAS_FORTRAN)
    begin_tests(netlib/odepack/test)
    set(odepacktests DLSODES DLSODAR DLSODI DLSODPK
      DLSODA DLSODE DLSODIS DLSODKR DLSOIBT)
    foreach(odetest IN LISTS odepacktests)
      new_test(SOURCES ${odetest}-test.f DEPS BLAS::BLAS)
      target_compile_options(${odetest}-test PRIVATE "-w")
      target_compile_options(${odetest}-test PRIVATE "-fallow-argument-mismatch")
      set_property(TARGET ${odetest}-test PROPERTY LINKER_LANGUAGE Fortran)

    endforeach()
    if(WITH_CXX)
      new_test(NAME odepacktest10 SOURCES test-funcC-inC.cpp# funC.cpp
        DEPS "LAPACK::LAPACK")
      target_compile_options(odepacktest10 PRIVATE "-w")
      #target_compile_features(odepacktest10 PUBLIC cxx_std_20)
      set_property(TARGET odepacktest10 PROPERTY LINKER_LANGUAGE CXX)
    endif(WITH_CXX)

    begin_tests(hairer/test)
    set(hairertests
      dr_iso
      dr_isosp
      #dr1_radau5
      #dr2_radau5
      #dr_radau
      #dr_radaup
      #dr_rodas
      # dr_seulex
      )
    foreach(hairertest IN LISTS hairertests)
      new_test(SOURCES ${hairertest}.f)
      target_compile_options(${hairertest} PRIVATE "-w")
      set_property(TARGET ${hairertest} PROPERTY LINKER_LANGUAGE Fortran)
    endforeach()
  
  if(WITH_MA57)
      begin_tests(lbl/example)
      new_test(NAME test_ma57 SOURCES example.c)
      configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${SOURCE_DIR}/lbl/example/mat.txt ${CMAKE_CURRENT_BINARY_DIR}/${SOURCE_DIR}/lbl/example/mat.txt COPYONLY)
    endif()
  endif()
  

endif()
