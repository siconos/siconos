# Makefile.am --- SICONOS

include $(top_srcdir)/config/Make-rules

check_PROGRAMS = test1
TESTS = $(check_PROGRAMS)
test1_SOURCES = SiconosMemoryTest.cpp TestMain.cpp
include_HEADERS = SiconosMemoryTest.h  
test1_INCLUDES = $(CPPUNIT_INCLUDES)
test1_LDFLAGS = $(CPPUNIT_LIBRARIES)
LIBS = $(SICONOSMATRIXVECTOR_LIB) $(SICONOSMEMORY_LIB) $(SICONOSXML_LIB) \
	$(SICONOSEXCEPTION_LIB) $(SICONOS_EXT_LDFLAGS) $(MODEL_LIB)
