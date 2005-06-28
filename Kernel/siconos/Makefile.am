# Makefile.am --- SICONOS

pkgdata_DATA = MakefilePluginSample

bin_SCRIPTS = siconos
siconos:siconos.in Makefile
	sed -e 's?\@prefix\@?$(prefix)?g' \
            -e 's/whatever//' \
	    < $(srcdir)/siconos.in > siconos.tmp
	sed -e 's?\@numericspath\@?$(NUMERICS_PATH)?g' \
            -e 's/whatever//' \
	    < $(srcdir)/siconos.tmp > siconos.tmp2 \
	&& mv siconos.tmp2 siconos
	chmod +x siconos

MakefilePluginSample:MakefilePluginSample.in Makefile
	sed -e 's?\@prefix\@?$(prefix)?g' \
            -e 's/whatever//' \
	    < $(srcdir)/MakefilePluginSample.in >MakefilePluginSample.tmp \
	&& mv MakefilePluginSample.tmp MakefilePluginSample
