# Makefile.am --- SICONOS

include $(top_srcdir)/config/Make-rules

SUBDIRS = src
EXTRA_DIST = 


doc:
	$(MKDIR) $(top_srcdir)/doc
	$(MKDIR) $(top_srcdir)/doc/api
	$(DOXYGEN) $(DOXYGEN_CFG)


docclean:
	$(RM) -rf $(top_srcdir)/doc/api
	$(RM) -rf $(top_srcdir)/doc

