#ifndef NONLINEARRELATION_CPP
#define NONLINEARRELATION_CPP

#include "NonlinearRelation.h"

//#include "const.h"


// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "debug.h"

NonlinearRelation::NonlinearRelation():
  FirstOrderType2R()
{
}

void NonlinearRelation::initializeWorkVectorsAndMatrices(Interaction& inter, VectorOfBlockVectors& DSlink, VectorOfVectors& workV, VectorOfSMatrices& workM)
{
  FirstOrderType2R::initializeWorkVectorsAndMatrices(inter, DSlink, workV, workM);
}

/*y = h(X)*/
void NonlinearRelation::computeh(double t, SiconosVector& x, SiconosVector& lambda, SiconosVector& y)
{
  DEBUG_PRINTF("NonlinearRelation::computeh at time %e\n ", t);
  DEBUG_EXPR(x.display());
  DEBUG_EXPR(lambda.display());
  y.setValue(0, 4.0 - x(0));
  y.setValue(1, 4.0 - x(1));
  y.setValue(2, 8.0 - x(0));
  y.setValue(3, 8.0 - x(1));
  DEBUG_EXPR(y.display());
}

/*g=g(lambda)*/
void NonlinearRelation::computeg(double t, SiconosVector& lambda, SiconosVector& r)
{
  DEBUG_PRINTF("NonlinearRelation::computeg at time %e\n ", t);
  DEBUG_EXPR(lambda.display());

  r(0) = 40.0 * (1 - lambda(2)) * (lambda(1));
  r(1) = 40.0 * (lambda(0)) * (1 - lambda(3));
  /*
  #ifdef SICONOS_DEBUG
    std::cout<<"NonlinearRelation::computeg with lambda="<<std::endl;
    lambda.display();
    std::cout<<std::endl;
    std::cout<<"NonlinearRelation::computeg modif g_alpha : \n";
    inter.data(g_alpha)->display();
    std::cout<<std::endl;
  #endif
  */
  DEBUG_EXPR(r.display());


}

void NonlinearRelation::computeJachlambda(double t, SiconosVector& x, SiconosVector& lambda, SimpleMatrix& D)
{
  DEBUG_PRINTF("NonlinearRelation::computeJachlambda at time %e\n ", t);
  D.zero();
}

void NonlinearRelation::computeJachx(double t, SiconosVector& x, SiconosVector& lambda, SimpleMatrix& C)
{
  DEBUG_PRINTF("NonlinearRelation::computeJachx at time %e\n ", t);

  C.setValue(0, 0, -1);
  C.setValue(0, 1, 0);
  C.setValue(1, 0, 0);
  C.setValue(1, 1, -1);
  C.setValue(2, 0, -1);
  C.setValue(2, 1, 0);
  C.setValue(3, 0, 0);
  C.setValue(3, 1, -1);
  DEBUG_EXPR(C.display());

}

void NonlinearRelation::computeJacglambda(double t, SiconosVector& lambda, SimpleMatrix& B)
{
  DEBUG_PRINTF("NonlinearRelation::computeJacglambda at time %e\n ", t);
  DEBUG_EXPR(lambda.display());


  B.setValue(0, 0, 0);
  B.setValue(1, 0, 40.0 * (1 - lambda(3)));

  B.setValue(0, 1, 40.0 * (1 - lambda(2)));
  B.setValue(1, 1, 0);

  B.setValue(0, 2, -40.0 * lambda(1));
  B.setValue(1, 2, 0);

  B.setValue(0, 3, 0);
  B.setValue(1, 3, -40.0 * lambda(0));
  DEBUG_EXPR(B.display());

}
#endif

