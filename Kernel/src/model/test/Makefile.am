# Makefile.am --- SICONOS

include $(top_srcdir)/config/Make-rules


check_PROGRAMS = test1
TESTS = $(check_PROGRAMS)
test1_SOURCES = PlatformTest.cpp TestMain.cpp
include_HEADERS = PlatformTest.h 
test1_INCLUDES = $(CPPUNIT_INCLUDES)
test1_LDFLAGS = $(CPPUNIT_LIBRARIES)
