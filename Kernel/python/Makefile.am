# Makefile.am --- SICONOS

include $(top_srcdir)/config/Make-rules

EXTRA_DIST = pySiconos.i createLibSiconos.sh

wrap: libsiconos runswig pySiconos_wrap.o lib

libsiconos:
	./createLibSiconos.sh
	g++ -g -shared -o libSiconos.so ./foo/*.o $(SICONOS_LDFLAGS) $(SICONOS_LDADD) $(LIBS)
	rm -rf ./foo

runswig: pySiconos.i
	swig -python -c++ pySiconos.i;

pySiconos_wrap.o: pySiconos_wrap.cxx
	./correctif.sh;
	g++ -c $(CXXFLAGS) $(INCLUDES) pySiconos_wrap.cxx -I/usr/include/python2.3;	

lib:
	g++ -shared $(SICONOS_LDFLAGS) $(SICONOS_LDADD) $(LIBS) libSiconos.so pySiconos_wrap.o -o _pySiconos.so 

clean:
	rm -f *.o *.so *.cxx *.py *.pyc
