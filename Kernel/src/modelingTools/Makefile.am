# Makefile.am --- SICONOS - ModelingTools

include $(top_srcdir)/config/Make-rules

.cpp.o:
	$(CXXCOMPILE) -c $<

SUBDIRS = test 

noinst_LTLIBRARIES = libModelingTools.la
libModelingTools_la_SOURCES =  BoundaryCondition.cpp ComplementarityConditionNSL.cpp \
		DSInputOutput.cpp DynamicalSystem.cpp EqualityConstraint.cpp Interaction.cpp \
		InteractionLink.cpp Topology.cpp \
		LagrangianDSIO.cpp LagrangianEC.cpp LagrangianLinearDSIO.cpp LagrangianLinearEC.cpp \
		LagrangianLinearR.cpp LagrangianDS.cpp LagrangianNonLinearR.cpp \
		LagrangianLinearTIDS.cpp LinearBC.cpp LinearDSIO.cpp LinearDS.cpp \
		LinearEC.cpp LinearTIEC.cpp LinearTIR.cpp NewtonImpactLawNSL.cpp NewtonImpactFrictionNSL.cpp NLinearBC.cpp \
		NonSmoothLaw.cpp NonSmoothDynamicalSystem.cpp PeriodicBC.cpp Relation.cpp RelayNSL.cpp \
		BoundaryCondition.h ComplementarityConditionNSL.h EqualityConstraint.h \
		DSInputOutput.h DynamicalSystem.h Interaction.h LagrangianDSIO.h LagrangianEC.h LagrangianLinearDSIO.h \
		LagrangianLinearEC.h LagrangianLinearR.h LagrangianDS.h LagrangianNonLinearR.h LagrangianLinearTIDS.h \
		LinearBC.h LinearDSIO.h LinearEC.h LinearDS.h LinearTIEC.h LinearTIR.h \
		NewtonImpactLawNSL.h NewtonImpactFrictionNSL.h NLinearBC.h NonSmoothLaw.h NonSmoothDynamicalSystem.h \
		PeriodicBC.h Relation.h RelayNSL.h
libModelingTools_la_LDFLAGS= -version-info $(SICONOS_SO_CURRENT):$(SICONOS_SO_REVISION):$(SICONOS_SO_AGE)

include_HEADERS = BoundaryCondition.h ComplementarityConditionNSL.h EqualityConstraint.h \
		DSInputOutput.h DynamicalSystem.h Interaction.h \
		InteractionLink.h Topology.h \
		LagrangianDSIO.h LagrangianEC.h LagrangianLinearDSIO.h LagrangianLinearEC.h LagrangianLinearR.h \
		LagrangianDS.h LagrangianNonLinearR.h LagrangianLinearTIDS.h LinearBC.h LinearDSIO.h \
		LinearEC.h LinearDS.h LinearTIEC.h LinearTIR.h NewtonImpactLawNSL.h NewtonImpactFrictionNSL.h \
		NLinearBC.h NonSmoothLaw.h NonSmoothDynamicalSystem.h PeriodicBC.h Relation.h RelayNSL.h
