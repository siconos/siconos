# Makefile.am --- SICONOS

include $(top_srcdir)/config/Make-rules
# include $(SICONOSPATH)/share/SICONOS/Make-rules
EXTRA_DIST = pySiconos.i createLibSiconos.sh

all:swig_interface lib

swig_interface:pySiconos.i

pySiconos_wrap.cxx:pySiconos.i
	$(SWIG) $(SWIG_PYTHON_OPT) -I$(SICONOS_INSTALLED_INCLUDES) pySiconos.i;
	./correctif.sh;

pySiconos_wrap.o:pySiconos_wrap.cxx
	g++ -c $(CXXFLAGS) -I$(SICONOS_INSTALLED_INCLUDES) $(SICONOS_EXT_INCLUDES) pySiconos_wrap.cxx $(PYTHON_CPPFLAGS)

#wrap: libsiconos runswig pySiconos_wrap.o lib

#libsiconos:
#	./createLibSiconos.sh
#	g++ -g -shared -o libSiconos.so ./foo/*.o $(SICONOS_LDFLAGS) $(SICONOS_LDADD) $(LIBS)
#	rm -rf ./foo

#runswig: pySiconos.i
#	swig -python -c++ pySiconos.i;

#pySiconos_wrap.o: pySiconos_wrap.cxx
#	./correctif.sh;
#	g++ -c $(CXXFLAGS) $(INCLUDES) pySiconos_wrap.cxx -I/usr/include/python2.3;	

lib:_pySiconos.so

_pySiconos.so:pySiconos_wrap.o
	$(CXX) $(CXX_DYNLIB)  $(SICONOS_INSTALLED_LIBS) $(SICONOS_EXT_LDFLAGS) pySiconos_wrap.o -o _pySiconos.so

#	g++ -shared $(SICONOS_LDFLAGS) $(SICONOS_LDADD) $(LIBS) libSiconos.so pySiconos_wrap.o -o _pySiconos.so 

clean:
	rm -f pySiconos_wrap.o _pySiconos.so pySiconos_wrap.cxx pySiconos.py pySiconos.pyc
