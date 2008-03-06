include $(top_srcdir)/config/Make-rules
TESTS_ENVIRONMENT = LD_LIBRARY_PATH=$(top_srcdir)/src/plugin/:$(top_builddir)/src/plugin/test:${LD_LIBRARY_PATH}

check_PROGRAMS = MainTest
TESTS = $(check_PROGRAMS)
MainTest_SOURCES = BouncingBall.cpp DiodeBridge.cpp BallBowl.cpp FunctionTest.cpp GlobalTests.cpp MainTest.cpp \
		    MainTest.h FunctionTest.h GlobalTests.h
siconos_HEADERS = MainTest.h FunctionTest.h GlobalTests.h
MainTest_INCLUDES = $(CPPUNIT_CFLAGS)
MainTest_LDFLAGS =  $(CPPUNIT_LIBS)  $(SICONOS_LIBS) $(SICONOS_EXT_LIBS)

pkglib_LTLIBRARIES = gTestPlugin.la
gTestPlugin_la_SOURCES = gTestPlugin.cpp
gTestPlugin_la_LDFLAGS = -module

all : gTestPlugin.so

.o.so:
	$(CXX) $(CXXFLAGS) $(CXX_DYNLIB) $< -o $@
