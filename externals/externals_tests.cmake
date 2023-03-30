include(tools4tests)

if(WITH_TESTING)

  add_custom_target(externals-tests echo "Start externals tests")
  
  if(HAS_FORTRAN)
    begin_tests(netlib/odepack/test) # DEPS "LAPACK::LAPACK")
    set(odepacktests DLSODES DLSODAR DLSODI DLSODPK
      DLSODA DLSODE DLSODIS DLSODKR DLSOIBT)
    foreach(odetest IN LISTS odepacktests)
      new_test(SOURCES ${odetest}-test.f)
      target_compile_options(${odetest}-test PRIVATE "-w")
      target_compile_options(${odetest}-test PRIVATE "-fallow-argument-mismatch")
    endforeach()
    if(WITH_CXX)
      new_test(NAME odepacktest10 SOURCES test-funcC-inC.cpp# funC.cpp
        DEPS "LAPACK::LAPACK")
      target_compile_options(odepacktest10 PRIVATE "-w")
      target_compile_features(odepacktest10 PUBLIC cxx_std_20)
    endif(WITH_CXX)
    
    begin_tests(hairer/test)# DEPS "LAPACK::LAPACK")
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
    #   new_test(NAME dr_iso1sp SOURCES dr_isosp.f)
    # target_compile_options(dr_iso1sp PRIVATE "-w")
    # new_test(SOURCES dr1_radau5.f)
    # new_test(SOURCES dr2_radau5.f)
    # target_compile_options(dr2_radau5 PRIVATE "-w")
    # new_test(SOURCES dr_radau.f)
    # new_test(SOURCES dr_radaup.f)
    # target_compile_options(dr_radaup PRIVATE "-w")
    # new_test(SOURCES dr_rodas.f)
    # target_compile_options(dr_rodas PRIVATE "-w")
    # new_test(SOURCES dr_seulex.f)
    # target_compile_options(dr_seulex PRIVATE "-w")
  endforeach()
  
  if(WITH_MA57)
      begin_tests(lbl/example)
      new_test(NAME test_ma57 SOURCES example.c)
      configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${SOURCE_DIR}/lbl/example/mat.txt ${CMAKE_CURRENT_BINARY_DIR}/${SOURCE_DIR}/lbl/example/mat.txt COPYONLY)
    endif()
  endif()
  

endif()
