# Makefile.am --- SICONOS

include $(top_srcdir)/config/Make-rules

SUBDIRS = config  ext src siconos sample #src must be compiled before main program
EXTRA_DIST = siconos.csh doc


documentation:
#	$(MKDIR) $(top_srcdir)/doc
	$(MKDIR) $(top_srcdir)/doc/api
	$(DOXYGEN) $(DOXYGEN_CFG)
   
documentationclean:
	$(RM) -rf $(top_srcdir)/doc/api

