# Makefile.am --- SICONOS

include $(top_srcdir)/config/Make-rules

.cpp.o:
	$(CXXCOMPILE) -c $<

lib_LIBRARIES = libModelStrategy.a
libModelStrategy_a_SOURCES = Adams.cpp EventDriven.cpp LCP.cpp Lsodar.cpp Moreau.cpp OneStepIntegrator.cpp OneStepNSProblem.cpp QP.cpp \
		Relay.cpp Strategy.cpp TimeDiscretisation.cpp TimeStepping.cpp \
		Adams.h EventDriven.h LCP.h Lsodar.h Moreau.h OneStepIntegrator.h OneStepNSProblem.h QP.h \
		Relay.h Strategy.h TimeDiscretisation.h TimeStepping.h	
include_HEADERS = Adams.h EventDriven.h LCP.h Lsodar.h Moreau.h OneStepIntegrator.h OneStepNSProblem.h QP.h \
		Relay.h Strategy.h TimeDiscretisation.h TimeStepping.h


