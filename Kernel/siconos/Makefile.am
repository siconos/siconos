# Makefile.am --- SICONOS

pkgdata_DATA = MakefilePluginSample

bin_SCRIPTS = siconos
siconos:siconos.in Makefile
	sed -e 's?\@prefix\@?$(prefix)?g' \
            -e 's/whatever//' \
	    < $(srcdir)/siconos.in > siconos.tmp \
	&& mv siconos.tmp siconos
	chmod +x siconos

MakefilePluginSample:Makefile
	sed -e 's?\@prefix\@?$(prefix)?g' \
            -e 's/whatever//' \
	    < $(srcdir)/MakefilePluginSample >MakefilePluginSample.tmp \
	&& mv MakefilePluginSample.tmp MakefilePluginSample
