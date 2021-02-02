include make.inc

# There is nothing to customize in this file. Customize make.inc.

# See if we are running on Mac.
ifneq ($(strip $(shell uname -a | grep -i 'darwin')),)
    PLATFORM = Mac
else
    PLATFORM = Linux
endif

SHELL = /bin/bash
TOP = ${PWD}
INCLUDEDIR = $(TOP)/include
CFLAGS += -I$(INCLUDEDIR)
LBL_LIBDIR  = $(TOP)/lib
LBL_OBJDIR  = $(TOP)/.objects
FCOMPILE = $(FC) $(FFLAGS)
CCOMPILE = $(CC) $(CFLAGS)
RM = rm -f
RMDIR = rmdir
SYMLINK = ln -sf

# On OSX, use Accelerate Framework if no other BLAS/LAPACK is specified.

ifeq ($(strip $(PLATFORM)),Mac)
ifeq ($(strip $(LIBBLAS)),)
ifeq ($(strip $(LIBLAPACK)),)
	LIBBLAS = -framework Accelerate
	LIBLAPACK =
endif
endif
endif

######################################################################

.SUFFIXES: .c .f

.PHONY: all clean veryclean mrclean distclean purge example

$(LBL_OBJDIR)/%.o: %.c
	$(CCOMPILE) -o $@ -c $<

$(LBL_OBJDIR)/%.o: %.f
	$(FCOMPILE) -o $@ -c $<

$(LBL_LIBDIR)/%.dylib: $(LBL_LIBDIR)/%.so
	$(SYMLINK) $^ $@

######################################################################

$(LBL_OBJDIR):
	[[ ! -d $(LBL_OBJDIR) ]] && mkdir $(LBL_OBJDIR)

$(LBL_LIBDIR):
	[[ ! -d $(LBL_LIBDIR) ]] && mkdir $(LBL_LIBDIR)

all:: $(LBL_OBJDIR) $(LBL_LIBDIR)

######################################################################

VPATH = $(HSL_DIR) $(TOP)/src

HSL_SRC = fd15d.f ma27ad.f mc21d.f mc34d.f mc47d.f mc59d.f mc64d.f mc71d.f ma57d.f
HSL_SRCS = $(addprefix $(HSL_DIR)/,$(HSL_SRC))
HSL_OBJS = $(addprefix $(LBL_OBJDIR)/,$(subst .f,.o,$(HSL_SRC)))

LBL_OBJS = lbl_mem.o ma27_lib.o ma57_lib.o lbl_lib.o
LBL_DEPS = $(addprefix $(LBL_OBJDIR)/,$(LBL_OBJS))

all:: liblbl

liblbl:: $(LBL_LIBDIR)/liblbl.so

$(LBL_LIBDIR)/liblbl.so:: $(LBL_DEPS) $(HSL_OBJS)
	$(FCOMPILE) $(NOFORMAIN) $(F_SHAR_LIB) -o $@ $+ $(LIBBLAS) $(LIBLAPACK) $(LIBMETIS) $(FLIBS)

ifeq ($(strip $(PLATFORM)),Mac)
liblbl:: $(LBL_LIBDIR)/liblbl.dylib
endif

######################################################################

example:: $(TOP)/example/example

$(TOP)/example/example: $(TOP)/example/example.c
	$(CCOMPILE) -o $@ $< -L$(TOP)/lib -llbl

######################################################################

clean:
	$(RM) $(LBL_OBJDIR)/*

veryclean: clean
	$(RM) $(LBL_LIBDIR)/*
	$(RMDIR) $(LBL_OBJDIR)
	$(RMDIR) $(LBL_LIBDIR)
	$(RM) $(TOP)/example/example

distclean: veryclean

mrclean: veryclean

purge : veryclean

######################################################################
