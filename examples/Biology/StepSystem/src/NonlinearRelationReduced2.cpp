#ifndef NONLINEARRELATIONREDUCED2_CPP
#define NONLINEARRELATIONREDUCED2_CPP

#include "NonlinearRelationReduced2.h"

//#include "const.h"
#define SICONOS_DEBUG

NonlinearRelationReduced2::NonlinearRelationReduced2():
  FirstOrderType2R()
{
}

void NonlinearRelationReduced2::initializeWorkVectorsAndMatrices(Interaction& inter, VectorOfBlockVectors& DSlink, VectorOfVectors& workV, VectorOfSMatrices& workM)
{
  FirstOrderType2R::initializeWorkVectorsAndMatrices(inter, DSlink, workV, workM);
}

/*y = h(X)*/
void NonlinearRelationReduced2::computeh(double t, SiconosVector& x, SiconosVector& lambda, SiconosVector& y)
{
  y(0) =  8.0 - x(0);
  y(1) = 8.0 - x(1);

}

/*g=g(lambda)*/
void NonlinearRelationReduced2::computeg(double t, SiconosVector& lambda, SiconosVector& r)
{
  r(0) = 40.0 * (1 - lambda(0));
  r(1) = 40.0 * (1 - lambda(1));
 }

void NonlinearRelationReduced2::computeJachlambda(double t, SiconosVector& x, SiconosVector& lambda, SimpleMatrix& D)
{
  D.zero();
}

void NonlinearRelationReduced2::computeJachx(double t, SiconosVector& x, SiconosVector& lambda, SimpleMatrix& C)
{

  C.setValue(0, 0, -1);
  C.setValue(0, 1, 0);
  C.setValue(1, 0, 0);
  C.setValue(1, 1, -1);

}


void NonlinearRelationReduced2::computeJacglambda(double t, SiconosVector& lambda, SimpleMatrix& B)
{

  B.setValue(0, 0, -40.0);
  B.setValue(1, 0,  0.0);
  B.setValue(0, 1,  0.0);
  B.setValue(1, 1, -40.0);
}
#endif

