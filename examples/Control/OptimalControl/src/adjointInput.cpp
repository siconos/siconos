#ifndef ADJOINTINPUT_CPP
#define ADJOINTINPUT_CPP

#include "adjointInput.hpp"
//#define SICONOS_DEBUG
adjointInput::adjointInput(): FirstOrderNonLinearR()
{
}


void adjointInput::initialize(Interaction& inter)
{
  FirstOrderNonLinearR::initialize(inter);
  K2.reset(new SimpleMatrix(2, 2));
  K2->setValue(0, 0, 0.0);
  K2->setValue(0, 1, -1.0 / 2.0);
  K2->setValue(1, 0, 1.0 / 2.0);
  K2->setValue(1, 1, 0.0);
}


double adjointInput::source(double t)
{
  double daux = 0;
  return daux;
}

/*y = h(X,lambda)*/
void adjointInput::computeh(double t, SiconosVector& x, SiconosVector& lambda, SiconosVector& z, SiconosVector& y)
{

#ifdef SICONOS_DEBUG
  std::cout << "********         computeH at " << t << std::endl;
#endif


  SP::SiconosVector betatmp(new SiconosVector(2));
  beta(t, x, betatmp);

  double betap = 2.0 * ((betatmp->getValue(0)) * x(2) + (betatmp->getValue(1)) * x(3));

  y(0) = lambda(1) + betap;
  y(1) =  2.0 - lambda(0);

#ifdef SICONOS_DEBUG
  std::cout << "modif heval : \n";
  y.display();
#endif

}



void adjointInput::computeg(double t, SiconosVector& x, SiconosVector& lambda, SiconosVector& z, SiconosVector& r)
{
#ifdef SICONOS_DEBUG
  std::cout << "************      computeG at: " << t << std::endl;
#endif

  SP::SiconosVector K2P(new SiconosVector(2));
  SP::SiconosVector P(new SiconosVector(2));
  P->setValue(0, x(2));
  P->setValue(1, x(3));



  prod(*K2, *P, *K2P, true);

  SP::SiconosVector betatmp(new SiconosVector(2));
  beta(t, x, betatmp);

  r(0) = betatmp->getValue(0) * (lambda(0) - 1.0);       //R=g_barre(x,lambda_barre)
  r(1) =  (betatmp->getValue(1)) * (lambda(0) - 1.0);
  r(2) = (K2P->getValue(0)) * (lambda(0) - 1.0);
  r(3) = (K2P->getValue(1)) * (lambda(0) - 1.0);


}

/** default function to compute jacobianH
 *  \param double : current time
 *  \param index for jacobian (0: jacobian according to x, 1 according to lambda)
 */
void adjointInput::computeJachx(double t, SiconosVector& x, SiconosVector& lambda, SiconosVector& z, SimpleMatrix& C)
{

#ifdef SICONOS_DEBUG
  std::cout << "computeJachx " << " at " << " " << t << std::endl;
#endif

  SP::SiconosVector betatmp(new SiconosVector(2));
  beta(t, x, betatmp);

  SP::SimpleMatrix jacbetaXtmp(new SimpleMatrix(2, 2));

  JacobianXbeta(t, x, jacbetaXtmp);


  C(0, 0) = x(3);
  C(0, 1) = -x(2) ;
  C(0, 2) = (-x(1) + 1.0) ;
  C(0, 3) = x(0);
  C(1, 0) = 0.0;
  C(1, 1) = 0.0 ;
  C(1, 2) = 0.0 ;
  C(1, 3) = 0.0;



#ifdef SICONOS_DEBUG
  std::cout << "modif Jachx : \n";
  C.display();
#endif


}
void adjointInput::computeJachlambda(double t, SiconosVector& x, SiconosVector& lambda, SiconosVector& z, SimpleMatrix& D)
{
#ifdef SICONOS_DEBUG
  std::cout << "computeJachlambda " << " at " << " " << t << std::endl;
#endif

  D(0, 0) = 0.0              ;
  D(0, 1) = 1.0 ;
  D(1, 0) = -1.0             ;
  D(1, 1) = 0.0 ;

#ifdef SICONOS_DEBUG
  std::cout << "modif Jachlambda : \n";
  D.display();
#endif

}
/** default function to compute jacobianG according to lambda
 *  \param double : current time
 *  \param index for jacobian: at the time only one possible jacobian => i = 0 is the default value .
 */
void adjointInput::computeJacgx(double t, SiconosVector& x, SiconosVector& lambda, SiconosVector& z, SimpleMatrix& K)
{

#ifdef SICONOS_DEBUG
  std::cout << "computeJacgx " << " at " << " " << t << std::endl;
#endif


  SP::SimpleMatrix jacbetaXtmp(new SimpleMatrix(2, 2));

  JacobianXbeta(t, x, jacbetaXtmp);

  K(0, 0) = jacbetaXtmp->getValue(0, 0) * (lambda(0) - 1.0) ;
  K(0, 1) = jacbetaXtmp->getValue(0, 1) * (lambda(0) - 1.0)  ;
  K(0, 2) = 0.0;
  K(0, 3) = 0.0;
  K(1, 0) = jacbetaXtmp->getValue(1, 0) * (lambda(0) - 1.0) ;
  K(1, 1) = jacbetaXtmp->getValue(1, 1) * (lambda(0) - 1.0)  ;
  K(1, 2) = 0.0  ;
  K(1, 3) = 0.0 ;
  K(2, 0) = 0.0;
  K(2, 1) = 0.0;
  K(2, 2) = K2->getValue(0, 0) * (lambda(0) - 1.0);
  K(2, 3) = K2->getValue(0, 1) * (lambda(0) - 1.0);
  K(3, 0) = 0.0;
  K(3, 1) = 0.0;
  K(3, 2) = K2->getValue(1, 0) * (lambda(0) - 1.0);
  K(3, 3) = K2->getValue(1, 1) * (lambda(0) - 1.0);

#ifdef SICONOS_DEBUG
  std::cout << "modif Jacgx : \n";
  K.display();
#endif


}
void adjointInput::computeJacglambda(double t, SiconosVector& x, SiconosVector& lambda, SiconosVector& z, SimpleMatrix& B)
{

  double *g = B.getArray();
#ifdef SICONOS_DEBUG
  std::cout << "computeJacglambda " << " at " << " " << t << std::endl;
#endif

  SP::SiconosVector K2P(new SiconosVector(2));
  SP::SiconosVector P(new SiconosVector(2));
  P->setValue(0, x(2));
  P->setValue(1, x(3));

  prod(*K2, *P, *K2P, true);

  SP::SiconosVector betatmp(new SiconosVector(2));
  beta(t, x, betatmp);

  g[0] = betatmp->getValue(0)  ;
  g[4] = 0.0;
  g[1] = betatmp->getValue(1)  ;
  g[5] = 0.0;
  g[2] = K2P->getValue(0)      ;
  g[6] = 0.0;
  g[3] = K2P->getValue(1)      ;
  g[7] = 0.0 ;



#ifdef SICONOS_DEBUG
  std::cout << "modif Jacglambda : \n";
  B.display();
#endif

}


void adjointInput::beta(double t, SiconosVector& xvalue, SP::SiconosVector beta)
{


  beta->setValue(0, -1.0 / 2.0 * xvalue(1) + 1.0 / 2.0) ;
  beta->setValue(1, 1.0 / 2.0 * xvalue(0)) ;
#ifdef SICONOS_DEBUG
  std::cout << "beta\n" << std::endl;;
  beta->display();
#endif

}

void adjointInput::JacobianXbeta(double t, SiconosVector& xvalue, SP::SimpleMatrix JacXbeta)
{


  JacXbeta->setValue(0, 0, 0.0) ;
  JacXbeta->setValue(0, 1, -1.0 / 2.0) ;
  JacXbeta->setValue(1, 0, 1.0 / 2.0) ;
  JacXbeta->setValue(1, 1, 0.0) ;
#ifdef SICONOS_DEBUG
  std::cout << "JacXbeta\n" << std::endl;;
  JacXbeta->display();
#endif

}


#endif
