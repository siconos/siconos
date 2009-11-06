DEBUG=-g
INC_EIGEN=-I/scratch/cadoux/eigen/2.0.1
OPT_SPARSE=-I/usr/include/suitesparse/ -DEIGEN_UMFPACK_SUPPORT
LNK_SPARSE=-lumfpack -lamd

all: test

test: test.o solverAC2.o
	g++ -o test $(LNK_SPARSE) test.o solverAC2.o

solverAC2.o: solverAC2.h solverAC2.cpp
	g++ -c $(DEBUG) $(INC_EIGEN) $(OPT_SPARSE) solverAC2.cpp

test.o: test.cpp solverAC2.h
	g++ -c $(DEBUG) $(INC_EIGEN) $(OPT_SPARSE) test.cpp

clean:
	rm *.o
	rm test
