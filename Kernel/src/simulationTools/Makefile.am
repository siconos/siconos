# Makefile.am --- SICONOS - simulationTools

include $(top_srcdir)/config/Make-rules

.cpp.o:
	$(CXXCOMPILE) -c $<


noinst_LTLIBRARIES = libSimulationTools.la
libSimulationTools_la_SOURCES = Adams.cpp CFD.cpp EventDriven.cpp LCP.cpp Lsodar.cpp Moreau.cpp OneStepIntegrator.cpp \
		OneStepNSProblem.cpp QP.cpp \
		Relay.cpp Strategy.cpp TimeDiscretisation.cpp TimeStepping.cpp \
		Adams.h CFD.h EventDriven.h LCP.h Lsodar.h Moreau.h OneStepIntegrator.h OneStepNSProblem.h QP.h \
		Relay.h Strategy.h TimeDiscretisation.h TimeStepping.h	
libSimulationTools_la_LDFLAGS= -version-info $(SICONOS_SO_CURRENT):$(SICONOS_SO_REVISION):$(SICONOS_SO_AGE)
include_HEADERS = Adams.h CFD.h EventDriven.h LCP.h Lsodar.h Moreau.h OneStepIntegrator.h OneStepNSProblem.h QP.h \
		Relay.h Strategy.h TimeDiscretisation.h TimeStepping.h


