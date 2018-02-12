#ifndef NONLINEARRELATIONWITHSIGNINVERSED_CPP
#define NONLINEARRELATIONWITHSIGNINVERSED_CPP

#include "NonlinearRelationWithSignInversed.hpp"

//#include "const.h"
//#define SICONOS_DEBUG

NonlinearRelationWithSignInversed::NonlinearRelationWithSignInversed():
  FirstOrderType2R()
{
}

/*y = h(X)*/
void NonlinearRelationWithSignInversed::computeh(double t, SiconosVector& x, SiconosVector& lambda, SiconosVector& y)
{


#ifdef SICONOS_DEBUG
  std::cout << "******** NonlinearRelationWithSignInversed::computeh computeh at " << t << std::endl;
#endif


  y.setValue(0, x(0) - 4);
  y.setValue(1, x(1) - 4);
  y.setValue(2, x(0) - 8);
  y.setValue(3, x(1) - 8);
#ifdef SICONOS_DEBUG
  std::cout << "modif heval : \n";
  y.display();
#endif
}

/*g=g(lambda)*/
void NonlinearRelationWithSignInversed::computeg(double t, SiconosVector& lambda, SiconosVector& r)
{

#ifdef SICONOS_DEBUG
  std::cout << "*** NonlinearRelationWithSignInversed::computeg     computeg at: " << t << std::endl;
#endif


  r.setValue(0, 10.0 * (1 + lambda(2)) * (1 - lambda(1)));
  r.setValue(1, 10.0 * (1 - lambda(0)) * (1 + lambda(3)));

#ifdef SICONOS_DEBUG
  std::cout << "NonlinearRelationWithSignInversed::computeg with lambda=" << std::endl;
  lambda.display();
  std::cout << std::endl;
  std::cout << "NonlinearRelationWithSignInversed::computeg modif g_alpha : \n";
  r.display();
  std::cout << std::endl;
#endif

}

void NonlinearRelationWithSignInversed::computeJachlambda(double t, SiconosVector& x, SiconosVector& lambda, SimpleMatrix& D)
{
  //    double *h = &(*_jachlambda)(0,0);
#ifdef SICONOS_DEBUG
  std::cout << "NonlinearRelationWithSignInversed::computeJachlambda " << " at " << " " << t << std::endl;
#endif
  D.zero();

}

void NonlinearRelationWithSignInversed::computeJachx(double t, SiconosVector& x, SiconosVector& lambda, SimpleMatrix& C)
{
  C.setValue(0, 0, 1);
  C.setValue(0, 1, 0);
  C.setValue(1, 0, 0);
  C.setValue(1, 1, 1);
  C.setValue(2, 0, 1);
  C.setValue(2, 1, 0);
  C.setValue(3, 0, 0);
  C.setValue(3, 1, 1);

#ifdef SICONOS_DEBUG
  std::cout << "NonlinearRelationWithSignInversed::computeJachx computeJachx " << " at " << " " << t << ":" << std::endl;
  C.display();
  std::cout << std::endl;
#endif

}

void NonlinearRelationWithSignInversed::computeJacglambda(double t, SiconosVector& lambda, SimpleMatrix& B)
{
  B.setValue(0, 0, 0);
  B.setValue(1, 0, -10.0 * (1 + lambda(3)));
  B.setValue(0, 1, -10.0 * (1 + lambda(2)));
  B.setValue(1, 1, 0);
  B.setValue(0, 2, 10.0 * (1 - lambda(1)));
  B.setValue(1, 2, 0);
  B.setValue(0, 3, 0);
  B.setValue(1, 3, 10.0 * (1 - lambda(0)));

#ifdef SICONOS_DEBUG
  std::cout << "NonlinearRelationWithSignInversed::computeJacgx " << " at " << " " << t << std::endl;
  B.display();
  std::cout << std::endl;
#endif

}
#endif

