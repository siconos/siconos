# Makefile.am --- SICONOS

include $(top_srcdir)/config/Make-rules

.cpp.o:
	$(CXXCOMPILE) -c $<

SUBDIRS = test
#EXTRA_DIST =

lib_LIBRARIES = libSiconosXml.a
libSiconosXml_a_SOURCES = AdamsXML.cpp BoundaryConditionXML.cpp \
		ComplementarityConditionNSLXML.cpp DSInputOutputXML.cpp DSXML.cpp LagrangianLinearRXML LagrangianLinearECXML.cpp\
		EqualityConstraintXML.cpp InteractionXML.cpp LagrangianDSIOXML.cpp LagrangianECXML.cpp\
		LagrangianLinearRXML.cpp LagrangianDSXML.cpp LagrangianNonLinearRXML.cpp \
		LagrangianLinearTIDSXML.cpp LCPXML.cpp LinearBCXML.cpp LinearDSIOXML.cpp \
		LinearSystemDSXML.cpp LinearECXML.cpp LinearTIECXML.cpp LinearTIRXML.cpp\
		LsodarXML.cpp MoreauXML.cpp NewtonImpactLawNSLXML.cpp NewtonImpactFrictionNSLXML.cpp NLinearBCXML.cpp \
		NonSmoothLawXML.cpp NSDSXML.cpp OneStepIntegratorXML.cpp OneStepNSProblemXML.cpp \
		PeriodicBCXML.cpp QPXML.cpp RelationXML.cpp RelayNSLXML.cpp SiconosDOMTreeTools.cpp \
		SiconosMemoryXML.cpp SiconosModelXML.cpp StrategyXML.cpp TimeDiscretisationXML.cpp \
		AdamsXML.h BoundaryConditionXML.h ComplementarityConditionNSLXML.h DSInputOutputXML.h \
		DSXML.h EqualityConstraintXML.h InteractionXML.h LagrangianDSIOXML.h LagrangianECXML.h LagrangianLinearECXML.h\
		LagrangianLinearRXML.h LagrangianDSXML.h LagrangianNonLinearRXML.h \
		LagrangianLinearTIDSXML.h LCPXML.h LinearBCXML.h LinearDSIOXML.h LinearECXML.h \
		LinearSystemDSXML.h LinearTIECXML.h LinearTIRXML.h LsodarXML.h MoreauXML.h \
		NewtonImpactLawNSLXML.h NewtonImpactFrictionNSLXML.h NLinearBCXML.h NonSmoothLawXML.h NSDSXML.h \
		OneStepIntegratorXML.h OneStepNSProblemXML.h PeriodicBCXML.h QPXML.h RelationXML.h \
		RelayNSLXML.h SiconosDOMTreeTools.h SiconosMemoryXML.h SiconosModelXML.h \
		StrategyXML.h TimeDiscretisationXML.h XMLTagsName.h

include_HEADERS = AdamsXML.h BoundaryConditionXML.h ComplementarityConditionNSLXML.h \
		DSInputOutputXML.h DSXML.h EqualityConstraintXML.h InteractionXML.h \
		LagrangianDSIOXML.h LagrangianECXML.h LagrangianLinearRXML.h LagrangianDSXML.h LagrangianLinearECXML.h \
		LagrangianNonLinearRXML.h LagrangianLinearTIDSXML.h LCPXML.h LinearBCXML.h \
		LinearDSIOXML.h LinearECXML.h LinearSystemDSXML.h LinearTIECXML.h LinearTIRXML.h \
		LsodarXML.h MoreauXML.h NewtonImpactLawNSLXML.h NewtonImpactFrictionNSLXML.h NLinearBCXML.h NonSmoothLawXML.h \
		NSDSXML.h OneStepIntegratorXML.h OneStepNSProblemXML.h PeriodicBCXML.h QPXML.h \
		RelationXML.h RelayNSLXML.h SiconosDOMTreeTools.h SiconosMemoryXML.h \
		SiconosModelXML.h StrategyXML.h TimeDiscretisationXML.h XMLTagsName.h
