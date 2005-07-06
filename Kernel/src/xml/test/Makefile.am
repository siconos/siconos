# Makefile.am --- SICONOS

include $(top_srcdir)/config/Make-rules

check_PROGRAMS = test1
TESTS = $(check_PROGRAMS)
test1_SOURCES = BoundaryConditionXMLTest.cpp LagrangianDSXMLTest.cpp OneStepIntegratorXMLTest.cpp SiconosMemoryXMLTest.cpp \
		DynamicalSystemXMLTest.cpp LinearDSXMLTest.cpp RelayNSLXMLTest.cpp InteractionXMLTest.cpp LinearTIRXMLTest.cpp \
		SiconosDOMTreeToolsTest.cpp TimeDiscretisationXMLTest.cpp TestMain.cpp
include_HEADERS = BoundaryConditionXMLTest.h LagrangianDSXMLTest.h OneStepIntegratorXMLTest.h SiconosMemoryXMLTest.h \
		DynamicalSystemXMLTest.h LinearDSXMLTest.h RelayNSLXMLTest.h InteractionXMLTest.h LinearTIRXMLTest.h \
		SiconosDOMTreeToolsTest.h TimeDiscretisationXMLTest.h
test1_INCLUDES = $(CPPUNIT_INCLUDES)
test1_LDFLAGS = $(CPPUNIT_LIBRARIES) $(SICONOS_LIBS) $(SICONOS_EXT_LDFLAGS)
