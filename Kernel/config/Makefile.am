# Makefile.am --- SICONOS

SUBDIRS = xmlschema
EXTRA_DIST = 

confdir=$(pkgdatadir)
conf_DATA = KernelConfig.sh
KernelConfig.sh:KernelConfig.sh.in Makefile
	sed -e 's?\@SICONOS_EXT_INCLUDES\@?$(SICONOS_EXT_INCLUDES)?g' \
	    -e 's?\@SICONOS_EXT_LDFLAGS\@?$(SICONOS_EXT_LDFLAGS)?g' \
	    -e 's?\@CXXFLAGS\@?$(CXXFLAGS)?g' \
            -e 's?\@CXX_DYNLIB\@?$(CXX_DYNLIB)?g' \
	       < $(srcdir)/KernelConfig.sh.in > KernelConfig.sh.tmp \
	&& mv KernelConfig.sh.tmp KernelConfig.sh
