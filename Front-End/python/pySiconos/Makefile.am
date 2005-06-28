# Makefile.am --- SICONOS


bin_SCRIPTS = pySiconos
pySiconos:pySiconos.in Makefile
	sed -e 's?\@kernelpath\@?$(KERNEL_PATH)?g' \
            -e 's/whatever//' \
	    < $(srcdir)/pySiconos.in > pySiconos.tmp
	sed -e 's?\@numericspath\@?$(NUMERICS_PATH)?g' \
            -e 's/whatever//' \
	    < $(srcdir)/pySiconos.tmp > pySiconos.tmp2 \
	&& mv pySiconos.tmp2 pySiconos
	chmod +x pySiconos

