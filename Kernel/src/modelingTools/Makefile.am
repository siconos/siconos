# Makefile.am --- SICONOS

include $(top_srcdir)/config/Make-rules

.cpp.o:
	$(CXXCOMPILE) -c $<

lib_LIBRARIES = libModelingTools.a
libModelingTools_a_SOURCES =  BoundaryCondition.cpp ComplementarityConditionNSL.cpp \
		DSInputOutput.cpp DynamicalSystem.cpp EqualityConstraint.cpp Interaction.cpp \
		LagrangianDSIO.cpp LagrangianEC.cpp LagrangianLinearEC.cpp LagrangianLinearR.cpp LagrangianDS.cpp LagrangianNonLinearR.cpp \
		LagrangianLinearTIDS.cpp LinearBC.cpp LinearDSIO.cpp LinearEC.cpp LinearSystemDS.cpp \
		LinearEC.cpp LinearTIEC.cpp LinearTIR.cpp NewtonImpactLawNSL.cpp NewtonImpactFrictionNSL.cpp NLinearBC.cpp \
		NonSmoothLaw.cpp NonSmoothDynamicalSystem.cpp PeriodicBC.cpp Relation.cpp RelayNSL.cpp \
		BoundaryCondition.h ComplementarityConditionNSL.h EqualityConstraint.h \
		DSInputOutput.h DynamicalSystem.h Interaction.h LagrangianDSIO.h LagrangianEC.h LagrangianLinearEC.h \
		LagrangianLinearR.h LagrangianDS.h LagrangianNonLinearR.h LagrangianLinearTIDS.h \
		LinearBC.h LinearDSIO.h LinearEC.h LinearSystemDS.h LinearTIEC.h LinearTIR.h \
		NewtonImpactLawNSL.h NewtonImpactFrictionNSL.h NLinearBC.h NonSmoothLaw.h NonSmoothDynamicalSystem.h \
		PeriodicBC.h Relation.h RelayNSL.h #NSDS.cpp NSDS.h

include_HEADERS = BoundaryCondition.h ComplementarityConditionNSL.h EqualityConstraint.h \
		DSInputOutput.h DynamicalSystem.h Interaction.h LagrangianDSIO.h LagrangianEC.h LagrangianLinearEC.h  LagrangianLinearR.h \
		LagrangianDS.h LagrangianNonLinearR.h LagrangianLinearTIDS.h LinearBC.h LinearDSIO.h \
		LinearEC.h LinearSystemDS.h LinearTIEC.h LinearTIR.h NewtonImpactLawNSL.h NewtonImpactFrictionNSL.cpp \
		NLinearBC.h NonSmoothLaw.h NonSmoothDynamicalSystem.h PeriodicBC.h Relation.h RelayNSL.h #NSDS.h
