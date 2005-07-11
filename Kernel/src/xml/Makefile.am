# Makefile.am --- SICONOS - xml

include $(top_srcdir)/config/Make-rules

.cpp.o:
	$(CXXCOMPILE) -c $<

#SUBDIRS = test


noinst_LTLIBRARIES = libSiconosXml.la
libSiconosXml_la_SOURCES = AdamsXML.cpp BoundaryConditionXML.cpp CFDXML.cpp \
		ComplementarityConditionNSLXML.cpp DSInputOutputXML.cpp DSXML.cpp \
		EqualityConstraintXML.cpp InteractionXML.cpp LagrangianDSIOXML.cpp LagrangianECXML.cpp LagrangianLinearDSIOXML.cpp \
		LagrangianLinearECXML.cpp LagrangianLinearRXML.cpp LagrangianDSXML.cpp LagrangianNonLinearRXML.cpp \
		LagrangianLinearTIDSXML.cpp LCPXML.cpp LinearBCXML.cpp LinearDSIOXML.cpp \
		LinearDSXML.cpp LinearECXML.cpp LinearTIECXML.cpp LinearTIRXML.cpp\
		LsodarXML.cpp MoreauXML.cpp NewtonImpactLawNSLXML.cpp NewtonImpactFrictionNSLXML.cpp NLinearBCXML.cpp \
		NonSmoothLawXML.cpp NSDSXML.cpp OneStepIntegratorXML.cpp OneStepNSProblemXML.cpp \
		PeriodicBCXML.cpp QPXML.cpp RelationXML.cpp RelayNSLXML.cpp SiconosDOMTreeTools.cpp \
		SiconosMemoryXML.cpp SiconosModelXML.cpp StrategyXML.cpp TimeDiscretisationXML.cpp \
		AdamsXML.h BoundaryConditionXML.h CFDXML.h ComplementarityConditionNSLXML.h DSInputOutputXML.h \
		DSXML.h EqualityConstraintXML.h InteractionXML.h LagrangianDSIOXML.h LagrangianECXML.h LagrangianLinearDSIOXML.h \
		LagrangianLinearECXML.h LagrangianLinearRXML.h LagrangianDSXML.h LagrangianNonLinearRXML.h \
		LagrangianLinearTIDSXML.h LCPXML.h LinearBCXML.h LinearDSIOXML.h LinearECXML.h \
		LinearDSXML.h LinearTIECXML.h LinearTIRXML.h LsodarXML.h MoreauXML.h \
		NewtonImpactLawNSLXML.h NewtonImpactFrictionNSLXML.h NLinearBCXML.h NonSmoothLawXML.h NSDSXML.h \
		OneStepIntegratorXML.h OneStepNSProblemXML.h PeriodicBCXML.h QPXML.h RelationXML.h \
		RelayNSLXML.h SiconosDOMTreeTools.h SiconosMemoryXML.h SiconosModelXML.h \
		StrategyXML.h TimeDiscretisationXML.h XMLTagsName.h

libSiconosXml_la_LDFLAGS= -version-info $(SICONOS_SO_CURRENT):$(SICONOS_SO_REVISION):$(SICONOS_SO_AGE)
include_HEADERS = AdamsXML.h BoundaryConditionXML.h CFDXML.h ComplementarityConditionNSLXML.h \
		DSInputOutputXML.h DSXML.h EqualityConstraintXML.h InteractionXML.h \
		LagrangianDSIOXML.h LagrangianECXML.h LagrangianLinearDSIOXML.h LagrangianLinearECXML.h \
		LagrangianLinearRXML.h LagrangianDSXML.h \
		LagrangianNonLinearRXML.h LagrangianLinearTIDSXML.h LCPXML.h LinearBCXML.h \
		LinearDSIOXML.h LinearECXML.h LinearDSXML.h LinearTIECXML.h LinearTIRXML.h \
		LsodarXML.h MoreauXML.h NewtonImpactLawNSLXML.h NewtonImpactFrictionNSLXML.h NLinearBCXML.h NonSmoothLawXML.h \
		NSDSXML.h OneStepIntegratorXML.h OneStepNSProblemXML.h PeriodicBCXML.h QPXML.h \
		RelationXML.h RelayNSLXML.h SiconosDOMTreeTools.h SiconosMemoryXML.h \
		SiconosModelXML.h StrategyXML.h TimeDiscretisationXML.h XMLTagsName.h
