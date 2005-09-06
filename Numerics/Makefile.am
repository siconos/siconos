# Makefile.am --- SICONOS

include $(top_srcdir)/config/Make-rules

# Test directories if cppunit is present.
if CPPUNIT_COND
       TEST_OPT = test
endif	

SUBDIRS = src $(TEST_OPT)


EXTRA_DIST = 


doc:
	$(MKDIR) $(top_srcdir)/doc
	$(MKDIR) $(top_srcdir)/doc/api
	$(DOXYGEN) $(DOXYGEN_CFG)


docclean:
	$(RM) -rf $(top_srcdir)/doc/api
	$(RM) -rf $(top_srcdir)/doc

