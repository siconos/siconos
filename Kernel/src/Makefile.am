# Makefile.am --- SICONOS - src

SUBDIRS =  utils xml plugin modelingTools simulationTools model plugin_test
EXTRA_DIST = 

lib_LTLIBRARIES = libSiconosKernel.la
libSiconosKernel_la_LDFLAGS= -version-info $(SICONOS_SO_CURRENT):$(SICONOS_SO_REVISION):$(SICONOS_SO_AGE) -no-undefined
libSiconosKernel_la_SOURCES=
libSiconosKernel_la_LIBADD = $(top_builddir)/src/model/libModel.la \
	$(top_builddir)/src/simulationTools/libSimulationTools.la \
	$(top_builddir)/src/modelingTools/libModelingTools.la \
	$(top_builddir)/src/xml/libSiconosXml.la \
	$(top_builddir)/src/utils/SiconosMemory/libSiconosMemory.la \
	$(top_builddir)/src/utils/SiconosAlgebra/libSiconosAlgebra.la  \
	$(top_builddir)/src/utils/SiconosSharedLibrary/libSiconosSharedLibrary.la \
	$(top_builddir)/src/utils/SiconosException/libSiconosException.la 
