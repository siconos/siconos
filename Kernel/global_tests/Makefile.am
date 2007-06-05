include $(top_srcdir)/config/Make-rules
TESTS_ENVIRONMENT = LD_LIBRARY_PATH=$(top_srcdir)/src/plugin/:$(top_builddir)/src/plugin/test:${LD_LIBRARY_PATH}

check_PROGRAMS = MainTest
TESTS = $(check_PROGRAMS)
MainTest_SOURCES = BouncingBall.cpp DiodeBridge.cpp BallBowl.cpp FunctionTest.cpp GlobalTests.cpp MainTest.cpp \
		    MainTest.h FunctionTest.h GlobalTests.h
siconos_HEADERS = MainTest.h FunctionTest.h GlobalTests.h
MainTest_INCLUDES = $(CPPUNIT_CFLAGS)
MainTest_LDFLAGS =  $(CPPUNIT_LIBS)  $(SICONOS_LIBS) $(SICONOS_EXT_LIBS)

all : gTestPlugin.so

gTestPlugin.o: gTestPlugin.cpp
	$(CXX) $(CXXFLAGS) -c gTestPlugin.cpp

gTestPlugin.so: gTestPlugin.o
	$(CXX) $(CXXFLAGS) $(CXX_DYNLIB) $< -o $@
