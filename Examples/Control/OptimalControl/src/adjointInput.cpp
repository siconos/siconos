#ifndef ADJOINTINPUT_CPP
#define ADJOINTINPUT_CPP

#include "adjointInput.hpp"
//#define SICONOS_DEBUG
adjointInput::adjointInput(): FirstOrderNonLinearR()
{
}


void adjointInput::initComponents(Interaction& inter, VectorOfBlockVectors& DSlink, VectorOfVectors& workV, VectorOfSMatrices& workM)
{
  FirstOrderNonLinearR::initComponents(inter, DSlink, workV, workM);

  SiconosVector& x = *workV[FirstOrderRVec::x];
  x = *DSlink[FirstOrderRDS::x];
  SiconosVector& lambda = *inter.lambda(0);

  K2.reset(new SimpleMatrix(2, 2));
  K2->setValue(0, 0, 0.0);
  K2->setValue(0, 1, -1.0 / 2.0);
  K2->setValue(1, 0, 1.0 / 2.0);
  K2->setValue(1, 1, 0.0);

  double t0 = 0;

  lambda.setValue(0, 0);
  lambda.setValue(1, 0);

  computeg(t0, x, lambda, *workV[FirstOrderRVec::g_alpha]);
  *DSlink[FirstOrderRDS::r] = *workV[FirstOrderRVec::g_alpha];
  computeJach(t0, inter, DSlink, workV, workM);
  computeJacg(t0, inter, DSlink, workV, workM);

}


double adjointInput::source(double t)
{
  double daux = 0;
  return daux;
}

/*y = h(X,lambda)*/
void adjointInput::computeh(double t, SiconosVector& x, SiconosVector& lambda, SiconosVector& y)
{

#ifdef SICONOS_DEBUG
  std::cout << "********         computeH at " << t << std::endl;
#endif


  SP::SiconosVector betatmp(new SiconosVector(2));
  beta(t, x, betatmp);

  double betap = 2.0 * ((betatmp->getValue(0)) * x(2) + (betatmp->getValue(1)) * x(3));



  y(0) = lambda(1) + betap;   //y_barre =Heval(x,lambda)
  y(1) =  2.0 - lambda(0);

#ifdef SICONOS_DEBUG
  std::cout << "modif heval : \n";
  y.display();
#endif

}



void adjointInput::computeg(double t, SiconosVector& x, SiconosVector& lambda, SiconosVector& r)
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
void adjointInput::computeJachx(double t, SiconosVector& x, SiconosVector& lambda, SimpleMatrix& C)
{

  double *h = C.getArray();
#ifdef SICONOS_DEBUG
  std::cout << "computeJachx " << " at " << " " << t << std::endl;
#endif

  SP::SiconosVector betatmp(new SiconosVector(2));
  beta(t, x, betatmp);

  SP::SiconosMatrix jacbetaXtmp(new SimpleMatrix(2, 2));

  JacobianXbeta(t, x, jacbetaXtmp);


  h[0] = x(3);
  h[2] = -x(2) ;
  h[4] = (-x(1) + 1.0) ;
  h[6] = x(0);
  h[1] = 0.0;
  h[3] = 0.0 ;
  h[5] = 0.0 ;
  h[7] = 0.0;



#ifdef SICONOS_DEBUG
  std::cout << "modif Jachx : \n";
  C.display();
#endif


}
void adjointInput::computeJachlambda(double t, SiconosVector& x, SiconosVector& lambda, SimpleMatrix& D)
{
  double *h = D.getArray();
#ifdef SICONOS_DEBUG
  std::cout << "computeJachlambda " << " at " << " " << t << std::endl;
#endif

  h[0] = 0.0              ;
  h[2] = 1.0 ;
  h[1] = -1.0             ;
  h[3] = 0.0 ;

#ifdef SICONOS_DEBUG
  std::cout << "modif Jachlambda : \n";
  D.display();
#endif

}
/** default function to compute jacobianG according to lambda
 *  \param double : current time
 *  \param index for jacobian: at the time only one possible jacobian => i = 0 is the default value .
 */
void adjointInput::computeJacgx(double t, SiconosVector& x, SiconosVector& lambda, SimpleMatrix& K)
{

  double *g = K.getArray();
#ifdef SICONOS_DEBUG
  std::cout << "computeJacgx " << " at " << " " << t << std::endl;
#endif


  SP::SiconosMatrix jacbetaXtmp(new SimpleMatrix(2, 2));

  JacobianXbeta(t, x, jacbetaXtmp);

  g[0] = jacbetaXtmp->getValue(0, 0) * (lambda(0) - 1.0) ;
  g[4] = jacbetaXtmp->getValue(0, 1) * (lambda(0) - 1.0)  ;
  g[8] = 0.0;
  g[12] = 0.0;
  g[1] = jacbetaXtmp->getValue(1, 0) * (lambda(0) - 1.0) ;
  g[5] = jacbetaXtmp->getValue(1, 1) * (lambda(0) - 1.0)  ;
  g[9] = 0.0  ;
  g[13] = 0.0 ;
  g[2] = 0.0;
  g[6] = 0.0;
  g[10] = K2->getValue(0, 0) * (lambda(0) - 1.0);
  g[14] = K2->getValue(0, 1) * (lambda(0) - 1.0);
  g[3] = 0.0;
  g[7] = 0.0;
  g[11] = K2->getValue(1, 0) * (lambda(0) - 1.0);
  g[15] = K2->getValue(1, 1) * (lambda(0) - 1.0);

#ifdef SICONOS_DEBUG
  std::cout << "modif Jacgx : \n";
  K.display();
#endif


}
void adjointInput::computeJacglambda(double t, SiconosVector& x, SiconosVector& lambda, SimpleMatrix& B)
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

void adjointInput::JacobianXbeta(double t, SiconosVector& xvalue, SP::SiconosMatrix JacXbeta)
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
