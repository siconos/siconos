#ifndef NONLINEARRELATION_CPP
#define NONLINEARRELATION_CPP

#include "NonlinearRelation.h"

//#include "const.h"
#define SICONOS_DEBUG

NonlinearRelation::NonlinearRelation():
  FirstOrderType2R()
{
}

void NonlinearRelation::initComponents(Interaction& inter, VectorOfBlockVectors& DSlink, VectorOfVectors& workV, VectorOfSMatrices& workM)
{
  FirstOrderType2R::initComponents(inter, DSlink, workV, workM);
}

/*y = h(X)*/
void NonlinearRelation::computeh(double t, SiconosVector& x, SiconosVector& lambda, SiconosVector& y)
{

  // SiconosVector& lambda = *inter.lambda(0);

  /*
  #ifdef SICONOS_DEBUG
    std::cout<<"******** NonlinearRelation::computeh computeh at "<<t<<std::endl;
  #endif
  */

  y.setValue(0, 4.0 - x(0));
  y.setValue(1, 4.0 - x(1));
  y.setValue(2, 8.0 - x(0));
  y.setValue(3, 8.0 - x(1));

}

/*g=g(lambda)*/
void NonlinearRelation::computeg(double t, SiconosVector& lambda, SiconosVector& r)
{
  /*
  #ifdef SICONOS_DEBUG
    std::cout<<"*** NonlinearRelation::computeg     computeg at: "<<t<<std::endl;
  #endif
  */

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

}

void NonlinearRelation::computeJachlambda(double t, SiconosVector& x, SiconosVector& lambda, SimpleMatrix& D)
{
  //    double *h = &(*_jachlambda)(0,0);
  /*
    #ifdef SICONOS_DEBUG
    std::cout<<"NonlinearRelation::computeJachlambda " <<" at " <<" "<<t<<std::endl;
    #endif
  */
  D.zero();
}

void NonlinearRelation::computeJachx(double t, SiconosVector& x, SiconosVector& lambda, SimpleMatrix& C)
{
  C.setValue(0, 0, -1);
  C.setValue(0, 1, 0);
  C.setValue(1, 0, 0);
  C.setValue(1, 1, -1);
  C.setValue(2, 0, -1);
  C.setValue(2, 1, 0);
  C.setValue(3, 0, 0);
  C.setValue(3, 1, -1);

}

void NonlinearRelation::computeJacglambda(double t, SiconosVector& lambda, SimpleMatrix& B)
{

  //  double *g = &(*Jacglambda)(0,0);
  B.setValue(0, 0, 0);
  B.setValue(1, 0, 40.0 * (1 - lambda(3)));

  B.setValue(0, 1, 40.0 * (1 - lambda(2)));
  B.setValue(1, 1, 0);

  B.setValue(0, 2, -40.0 * lambda(1));
  B.setValue(1, 2, 0);

  B.setValue(0, 3, 0);
  B.setValue(1, 3, -40.0 * lambda(0));

}
#endif

