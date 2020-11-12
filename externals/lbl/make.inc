# Specify location of MA27, MA57 and METIS.

HSL_DIR  = ext
LIBMETIS  = -L${HOME}/lib -lmetis

# On OSX, leave LIBBLAS and LIBLAPACK blank to use Accelerate.

LIBBLAS   = #-lblas
LIBLAPACK = #-llapack

# Specify compiler and options.

FC = gfortran
FFLAGS = -g #-O2
NOFORMAIN = #-nofor_main
F_SHAR_LIB = -shared
FLIBS =
CC = gcc
CFLAGS = -g #-O2
CLIBS =
