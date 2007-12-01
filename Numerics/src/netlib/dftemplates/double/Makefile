FC = f77

FFLAGS= -g -u

libs = dlapack.o dblas.o

tester_o  = Tester.o MatGen.o dspdchk.o dnsychk.o \
            Utils.o \
            MatVec.o PSolve.o \
			CG.o CGREVCOM.o \
			BiCG.o BiCGREVCOM.o \
			BiCGSTAB.o BiCGSTABREVCOM.o \
			CGS.o CGSREVCOM.o \
			Cheby.o ChebyREVCOM.o \
			GMRES.o GMRESREVCOM.o \
			Jacobi.o JacobiREVCOM.o \
			QMR.o QMRREVCOM.o \
			SOR.o SORREVCOM.o \
			STOPTEST2.o \
            $(libs)

.f.o:;  $(FC) $(FFLAGS) -c $<

xdtest:   $(tester_o)
	$(FC) $(FFLAGS) -o $@ $(tester_o)

clean:
	rm -f core a.out xdtest* *.trace *.o
