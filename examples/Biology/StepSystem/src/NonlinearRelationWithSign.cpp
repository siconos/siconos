#ifndef NONLINEARRELATIONWITHSIGN_CPP
#define NONLINEARRELATIONWITHSIGN_CPP

#include "NonlinearRelationWithSign.h"

//#include "const.h"
#define SICONOS_DEBUG

NonlinearRelationWithSign::NonlinearRelationWithSign():
  FirstOrderType2R()
{
}

void NonlinearRelationWithSign::initializeWorkVectorsAndMatrices(Interaction& inter, VectorOfBlockVectors& DSlink, VectorOfVectors& workV, VectorOfSMatrices& workM)
{
  FirstOrderType2R::initializeWorkVectorsAndMatrices(inter, DSlink, workV, workM);
}

/*y = h(X)*/
void NonlinearRelationWithSign::computeh(double t, SiconosVector& x, SiconosVector& lambda, SiconosVector& y)
{

  //SiconosVector& lambda = *inter.lambda(0);

#ifdef SICONOS_DEBUG
  std::cout << "******** NonlinearRelationWithSign::computeh computeh at " << t << std::endl;
#endif


  y.setValue(0, 4.0 - x(0));
  y.setValue(1, 4.0 - x(1));
  y.setValue(2, 8.0 - x(0));
  y.setValue(3, 8.0 - x(1));
#ifdef SICONOS_DEBUG
  std::cout << "modif heval : \n";
  y.display();
#endif
}


/*g=g(lambda)*/
void NonlinearRelationWithSign::computeg(double t, SiconosVector& lambda, SiconosVector& r)
{

#ifdef SICONOS_DEBUG
  std::cout << "*** NonlinearRelationWithSign::computeg     computeg at: " << t << std::endl;
#endif


  r.setValue(0, 10.0 * (1 - lambda(2)) * (1 + lambda(1)));
  r.setValue(1, 10.0 * (1 + lambda(0)) * (1 - lambda(3)));

#ifdef SICONOS_DEBUG
  std::cout << "NonlinearRelationWithSign::computeg with lambda=" << std::endl;
  lambda.display();
  std::cout << std::endl;
  std::cout << "NonlinearRelationWithSign::computeg modif g_alpha : \n";
  r.display();
  std::cout << std::endl;
#endif

}

void NonlinearRelationWithSign::computeJachx(double t, SiconosVector& x, SiconosVector& lambda, SimpleMatrix& C)
{
  C.setValue(0, 0, -1);
  C.setValue(0, 1, 0);
  C.setValue(1, 0, 0);
  C.setValue(1, 1, -1);
  C.setValue(2, 0, -1);
  C.setValue(2, 1, 0);
  C.setValue(3, 0, 0);
  C.setValue(3, 1, -1);

#ifdef SICONOS_DEBUG
  std::cout << "NonlinearRelationWithSign::computeJachx computeJachx " << " at " << " " << t << ":" << std::endl;
  C.display();
  std::cout << std::endl;
#endif

}

void NonlinearRelationWithSign::computeJachlambda(double t, SiconosVector& x, SiconosVector& lambda, SimpleMatrix& D)
{
  D.zero();
}

void NonlinearRelationWithSign::computeJacglambda(double t, SiconosVector& lambda, SimpleMatrix& B)
{

  //  double *g = &(*Jacglambda)(0,0);
  B.setValue(0, 0, 0);
  B.setValue(1, 0, 10.0 * (1 - lambda(3)));
  B.setValue(0, 1, 10.0 * (1 - lambda(2)));
  B.setValue(1, 1, 0);
  B.setValue(0, 2, -10.0 * (1 + lambda(1)));
  B.setValue(1, 2, 0);
  B.setValue(0, 3, 0);
  B.setValue(1, 3, -10.0 * (1 + lambda(0)));

#ifdef SICONOS_DEBUG
  std::cout << "NonlinearRelationWithSign::computeJacgx " << " at " << " " << t << std::endl;
  B.display();
  std::cout << std::endl;
#endif

}
#endif

