# Makefile.am --- SICONOS

include $(top_srcdir)/config/Make-rules

check_PROGRAMS = test1
TESTS = $(check_PROGRAMS)
test1_SOURCES = BoundaryConditionXMLTest.cpp LagrangianDSXMLTest.cpp OneStepIntegratorXMLTest.cpp SiconosMemoryXMLTest.cpp \
		DynamicalSystemXMLTest.cpp LinearSystemDSXMLTest.cpp RelayNSLXMLTest.cpp InteractionXMLTest.cpp LinearTIRXMLTest.cpp \
		SiconosDOMTreeToolsTest.cpp TimeDiscretisationXMLTest.cpp TestMain.cpp
include_HEADERS = BoundaryConditionXMLTest.h LagrangianDSXMLTest.h OneStepIntegratorXMLTest.h SiconosMemoryXMLTest.h \
		DynamicalSystemXMLTest.h LinearSystemDSXMLTest.h RelayNSLXMLTest.h InteractionXMLTest.h LinearTIRXMLTest.h \
		SiconosDOMTreeToolsTest.h TimeDiscretisationXMLTest.h  


test1_INCLUDES = $(CPPUNIT_INCLUDES)
test1_LDFLAGS = $(CPPUNIT_LIBRARIES) 

LIBS=$(MODEL_LIB) $(SIMULATIONTOOLS_LIB) $(MODELINGTOOLS_LIB) \
	$(SICONOSXML_LIB) $(SICONOSMEMORY_LIB) $(SICONOSMATRIXVECTOR_LIB) \
	$(SICONOSSHAREDLIBRARY_LIB) $(SICONOSEXCEPTION_LIB) $(SICONOS_EXT_LDFLAGS) \
	$(MODEL_LIB) $(SIMULATIONTOOLS_LIB) $(SICONOSNUMERICS_LIBRARIES)  \
	$(MODELINGTOOLS_LIB) $(SICONOSMEMORY_LIB) $(SICONOSALGEBRA_LIB) $(SICONOSSHAREDLIBRARY_LIB)
