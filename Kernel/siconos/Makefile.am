# Makefile.am --- SICONOS

pkgdata_DATA = MakefilePluginSample

bin_SCRIPTS = siconos
siconos:siconos.in
	sed -e 's/whatever//' < $(srcdir)/siconos.in > siconos
	chmod +x siconos

