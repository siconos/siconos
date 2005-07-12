# Makefile.am --- SICONOS

include $(top_srcdir)/config/Make-rules


check_PROGRAMS = test1
TESTS = $(check_PROGRAMS)
test1_SOURCES = ModelTest.cpp TestMain.cpp
include_HEADERS = ModelTest.h 
test1_INCLUDES = $(CPPUNIT_INCLUDES)
test1_LDFLAGS = $(CPPUNIT_LIBRARIES) $(SICONOS_LIBS)  $(SICONOS_EXT_LDFLAGS)
