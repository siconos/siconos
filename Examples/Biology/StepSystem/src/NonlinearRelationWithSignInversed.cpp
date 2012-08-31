#ifndef NONLINEARRELATIONWITHSIGNINVERSED_CPP
#define NONLINEARRELATIONWITHSIGNINVERSED_CPP

#include "NonlinearrelationWithSignInversed.h"

//#include "const.h"
#define SICONOS_DEBUG

NonlinearRelationWithSignInversed::NonlinearRelationWithSignInversed():
  FirstOrderType2R()
{
}

void NonlinearRelationWithSignInversed::initialize(Interaction& inter)
{
  FirstOrderType2R::initialize(inter);
  unsigned int sizeY = inter.getSizeOfY();
  unsigned int sizeDS = inter.getSizeOfDS();
  SiconosVector& y = *inter.y(0);
  SiconosVector& lambda = *inter.lambda(0);

  double t0 = 0;

  _jachx->resize(sizeY, sizeDS);
  _jachlambda->resize(sizeY, sizeY);

  _jacgx->resize(sizeDS, sizeDS);
  _jacglambda->resize(sizeDS, sizeY);


  //  SiconosVector workX = *inter.data(x);
  //_workX->setValue(0,0);
  y.setValue(0, 0);
  y.setValue(1, 0);
  y.setValue(2, 0);
  y.setValue(3, 0);
  lambda.setValue(0, 0);
  lambda.setValue(1, 0);
  lambda.setValue(2, 0);
  lambda.setValue(3, 0);


  //  computeH(t0);
  computeg(t0, inter);
  computeJach(t0, inter);
  computeJacg(t0, inter);
  *inter.data(r) = *inter.data(g_alpha);


#ifdef SICONOS_DEBUG
  std::cout << "data[r (g_alpha)] init\n";
  inter.data(r)->display();
#endif


}

/*y = h(X)*/
void NonlinearRelationWithSignInversed::computeh(double t, Interaction& inter)
{

  SiconosVector workX = *inter.data(x);
  //SiconosVector& lambda = *inter.lambda(0);

#ifdef SICONOS_DEBUG
  std::cout << "******** NonlinearRelationWithSignInversed::computeh computeh at " << t << std::endl;
#endif

  SP::SiconosVector Heval = inter.Halpha();

  Heval->setValue(0, workX(0) - 4);
  Heval->setValue(1, workX(1) - 4);
  Heval->setValue(2, workX(0) - 8);
  Heval->setValue(3, workX(1) - 8);
#ifdef SICONOS_DEBUG
  std::cout << "modif heval : \n";
  Heval->display();
#endif
}

/*g=g(lambda)*/
void NonlinearRelationWithSignInversed::computeg(double t, Interaction& inter)
{
  SiconosVector& lambda = *inter.lambda(0);

#ifdef SICONOS_DEBUG
  std::cout << "*** NonlinearRelationWithSignInversed::computeg     computeg at: " << t << std::endl;
#endif


  inter.data(g_alpha)->setValue(0, 10.0 * (1 + lambda(2)) * (1 - lambda(1)));
  inter.data(g_alpha)->setValue(1, 10.0 * (1 - lambda(0)) * (1 + lambda(3)));

#ifdef SICONOS_DEBUG
  std::cout << "NonlinearRelationWithSignInversed::computeg with lambda=" << std::endl;
  lambda.display();
  std::cout << std::endl;
  std::cout << "NonlinearRelationWithSignInversed::computeg modif g_alpha : \n";
  inter.data(g_alpha)->display();
  std::cout << std::endl;
#endif

}

void NonlinearRelationWithSignInversed::computeJachlambda(double t, Interaction& inter)
{
  //    double *h = &(*_jachlambda)(0,0);
#ifdef SICONOS_DEBUG
  std::cout << "NonlinearRelationWithSignInversed::computeJachlambda " << " at " << " " << t << std::endl;
#endif
  _jachlambda->zero();

}

void NonlinearRelationWithSignInversed::computeJachx(double t, Interaction& inter)
{
  _jachx->setValue(0, 0, 1);
  _jachx->setValue(0, 1, 0);
  _jachx->setValue(1, 0, 0);
  _jachx->setValue(1, 1, 1);
  _jachx->setValue(2, 0, 1);
  _jachx->setValue(2, 1, 0);
  _jachx->setValue(3, 0, 0);
  _jachx->setValue(3, 1, 1);

#ifdef SICONOS_DEBUG
  std::cout << "NonlinearRelationWithSignInversed::computeJachx computeJachx " << " at " << " " << t << ":" << std::endl;
  _jachx->display();
  std::cout << std::endl;
#endif

}

void NonlinearRelationWithSignInversed::computeJacgx(double t, Interaction& inter)
{
  //   double *g = &(*jacgx)(0,0);
  _jacgx->zero();

#ifdef SICONOS_DEBUG
  std::cout << "NonlinearRelationWithSignInversed::computeJacgx " << " at " << " " << t << std::endl;
  _jacgx->display();
  std::cout << std::endl;
#endif

}

void NonlinearRelationWithSignInversed::computeJacglambda(double t, Interaction& inter)
{

  SiconosVector& lambda = *inter.lambda(0);

  //  double *g = &(*Jacglambda)(0,0);
  _jacglambda->setValue(0, 0, 0);
  _jacglambda->setValue(1, 0, -10.0 * (1 + lambda(3)));
  _jacglambda->setValue(0, 1, -10.0 * (1 + lambda(2)));
  _jacglambda->setValue(1, 1, 0);
  _jacglambda->setValue(0, 2, 10.0 * (1 - lambda(1)));
  _jacglambda->setValue(1, 2, 0);
  _jacglambda->setValue(0, 3, 0);
  _jacglambda->setValue(1, 3, 10.0 * (1 - lambda(0)));

#ifdef SICONOS_DEBUG
  std::cout << "NonlinearRelationWithSignInversed::computeJacgx " << " at " << " " << t << std::endl;
  _jacglambda->display();
  std::cout << std::endl;
#endif

}
#endif

