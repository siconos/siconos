# Makefile.am --- SICONOS

include $(top_srcdir)/config/Make-rules

SUBDIRS = config src siconos #src must be compiled before main program
EXTRA_DIST = kernel.csh.TEMPLATE


documentation:
	$(MKDIR) $(top_srcdir)/doc
	$(MKDIR) $(top_srcdir)/doc/api
	$(DOXYGEN) $(DOXYGEN_CFG)
   
documentationclean:
	$(RM) -rf $(top_srcdir)/doc/api
	$(RM) -rf $(top_srcdir)/doc

