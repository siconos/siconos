# Makefile.am --- SICONOS

include $(top_srcdir)/config/Make-rules

#SUBDIRS = 
#EXTRA_DIST =

all : BasicPlugin.so

BasicPlugin.o: BasicPlugin.c
	$(CXX) $(CXXFLAGS) -c BasicPlugin.c

BasicPlugin.so: BasicPlugin.o
	$(CXX) $(CXXFLAGS) $(CXX_DYNLIB) $< -o $@

clean:
	$(RM) -f BasicPlugin.o BasicPlugin.so

pkgdata_DATA = BasicPlugin.so


