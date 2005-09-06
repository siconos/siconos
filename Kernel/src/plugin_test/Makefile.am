# Makefile.am --- SICONOS

include $(top_srcdir)/config/Make-rules

#SUBDIRS = 
#EXTRA_DIST =

all : TestPlugin.so

TestPlugin.o: TestPlugin.c
	$(CXX) $(CXXFLAGS) -c TestPlugin.c

TestPlugin.so: TestPlugin.o
	$(CXX) $(CXXFLAGS) $(CXX_DYNLIB) $< -o $@

clean:
	$(RM) -f TestPlugin.o TestPlugin.so

pkgdata_DATA = TestPlugin.so


