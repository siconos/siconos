# Makefile.am --- SICONOS

include $(top_srcdir)/config/Make-rules

.cpp.o:
	$(CXXCOMPILE) -c $<

lib_LIBRARIES = libModelFormalisation.a
libModelFormalisation_a_SOURCES =  BoundaryCondition.cpp ComplementarityConditionNSL.cpp \
		DSInputOutput.cpp DynamicalSystem.cpp EqualityConstraint.cpp Interaction.cpp \
		LagrangianDSIO.cpp LagrangianEC.cpp LagrangianLinearR.cpp LagrangianNLDS.cpp LagrangianNonLinearR.cpp \
		LagrangianTIDS.cpp LinearBC.cpp LinearDSIO.cpp LinearEC.cpp LinearSystemDS.cpp \
		LinearEC.cpp LinearTIEC.cpp LinearTIR.cpp NewtonImpactLawNSL.cpp NewtonImpactFrictionNSL.cpp NLinearBC.cpp \
		NonSmoothLaw.cpp NonSmoothDynamicalSystem.cpp PeriodicBC.cpp Relation.cpp RelayNSL.cpp \
		BoundaryCondition.h ComplementarityConditionNSL.h EqualityConstraint.h \
		DSInputOutput.h DynamicalSystem.h Interaction.h LagrangianDSIO.h LagrangianEC.h \
		LagrangianLinearR.h LagrangianNLDS.h LagrangianNonLinearR.h LagrangianTIDS.h \
		LinearBC.h LinearDSIO.h LinearEC.h LinearSystemDS.h LinearTIEC.h LinearTIR.h \
		NewtonImpactLawNSL.h NewtonImpactFrictionNSL.h NLinearBC.h NonSmoothLaw.h NonSmoothDynamicalSystem.h \
		PeriodicBC.h Relation.h RelayNSL.h #NSDS.cpp NSDS.h

include_HEADERS = BoundaryCondition.h ComplementarityConditionNSL.h EqualityConstraint.h \
		DSInputOutput.h DynamicalSystem.h Interaction.h LagrangianDSIO.h LagrangianEC.h LagrangianLinearR.h \
		LagrangianNLDS.h LagrangianNonLinearR.h LagrangianTIDS.h LinearBC.h LinearDSIO.h \
		LinearEC.h LinearSystemDS.h LinearTIEC.h LinearTIR.h NewtonImpactLawNSL.h NewtonImpactFrictionNSL.cpp \
		NLinearBC.h NonSmoothLaw.h NonSmoothDynamicalSystem.h PeriodicBC.h Relation.h RelayNSL.h #NSDS.h
