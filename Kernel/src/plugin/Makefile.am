# Makefile.am --- SICONOS

include $(top_srcdir)/config/Make-rules

#SUBDIRS = 
#EXTRA_DIST =

all : BasicPlugin.o

BasicPlugin.o: BasicPlugin.c
	$(CXX) $(CXXFLAGS) -c BasicPlugin.c
	$(CXX) $(CXXFLAGS) $(CXX_DYNLIB) -W -O BasicPlugin.o -o BasicPlugin.so

.c.so: 
	$(CXX) $(CXXFLAGS) -c $<
	$(CXX) $(CXXFLAGS) $(CXX_DYNLIB) *.o -o $@






