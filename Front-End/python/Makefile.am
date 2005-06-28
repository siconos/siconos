# Makefile.am --- SICONOS

SUBDIRS = pySiconos

include $(top_srcdir)/config/Make-rules
include $(KERNEL_PATH)/share/SICONOS/KernelConfig.sh
EXTRA_DIST = pySiconos.i createLibSiconos.sh

all:swig_interface lib

swig_interface:pySiconos.i

pySiconos_wrap.cxx:pySiconos.i
	$(SWIG) $(SWIG_PYTHON_OPT) -I$(SICONOS_INSTALLED_INCLUDES) pySiconos.i;
	sed "s/SWIG_append_msg/SWIG_append_errmsg/g" pySiconos_wrap.cxx > pySiconos_wrap.cxx.res ;
	mv pySiconos_wrap.cxx.res pySiconos_wrap.cxx;

pySiconos_wrap.o:pySiconos_wrap.cxx
	g++ -c $(CXXFLAGS) -I$(SICONOS_INSTALLED_INCLUDES) $(LIBEXT_INCLUDES) pySiconos_wrap.cxx $(PYTHON_CPPFLAGS)

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
	$(CXX) $(CXX_DYNLIB)  $(SICONOS_INSTALLED_LIBS) $(LIBEXT_LDFLAGS) pySiconos_wrap.o -o _pySiconos.so

#	g++ -shared $(SICONOS_LDFLAGS) $(SICONOS_LDADD) $(LIBS) libSiconos.so pySiconos_wrap.o -o _pySiconos.so 

clean:
	rm -f pySiconos_wrap.o _pySiconos.so pySiconos_wrap.cxx pySiconos.py pySiconos.pyc

install-data-am:
	@INSTALL@ -m 0644 _pySiconos.so $(PYTHON_SITE_PKG)
	@INSTALL@ -m 0644 pySiconos.py $(PYTHON_SITE_PKG)
