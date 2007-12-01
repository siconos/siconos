FC = gcc

FFLAGS =

f2cLIBS = F2CLIBS/libF77.a F2CLIBS/libI77.a

tester_o  = Tester.o MatGen.o dspdchk.o dnsychk.o \
            CG.o BiCG.o CGS.o BiCGSTAB.o GMRES.o QMR.o Cheby.o  \
            Jacobi.o SOR.o Utils.o \
            MatVec.o PSolve.o \
            dlapack.o dblas.o

.f.o:;  $(FC) $(FFLAGS) -c $<

xtest:   $(tester_o)
	$(FC) $(FFLAGS) -o $@ $(tester_o) $(f2cLIBS) -lm

clean:
	rm -f core a.out xtest* *.trace *.o
