#ifndef ADJOINTINPUT_CPP
#define ADJOINTINPUT_CPP

#include "adjointInput.hpp"
//#define SICONOS_DEBUG
adjointInput::adjointInput():
  FirstOrderType2R()
{
}


void adjointInput::initialize(Interaction& inter)
{
  FirstOrderType2R::initialize(inter);
  unsigned int sizeY = inter.getSizeOfY();
  unsigned int sizeDS = inter.getSizeOfDS();
  SiconosVector& lambda = *inter.lambda(0);

  K2 = new SimpleMatrix(2, 2);
  K2->setValue(0, 0, 0.0);
  K2->setValue(0, 1, -1.0 / 2.0);
  K2->setValue(1, 0, 1.0 / 2.0);
  K2->setValue(1, 1, 0.0);

  double t0 = 0;

  _jachx->resize(sizeY, sizeDS);
  _jachlambda->resize(sizeY, sizeY);

  _jacgx->resize(sizeDS, sizeDS);
  _jacglambda->resize(sizeDS, sizeY);


  lambda.setValue(0, 0);
  lambda.setValue(1, 0);

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


double adjointInput::source(double t)
{
  double daux = 0;
  return daux;
}

/*y = h(X,lambda)*/
void adjointInput::computeh(double t, Interaction& inter)
{

  SiconosVector workX = *inter.data(x);
  SiconosVector& lambda = *inter.lambda(0);

#ifdef SICONOS_DEBUG
  std::cout << "********         computeH at " << t << std::endl;
#endif
  SP::SiconosVector Heval = inter.Halpha();


  SP::SiconosVector betatmp(new SiconosVector(2));
  beta(t, workX, betatmp);

  double betap = 2.0 * ((betatmp->getValue(0)) * workX(2) + (betatmp->getValue(1)) * workX(3));



  Heval->setValue(0, lambda(1) + betap);   //y_barre =Heval(x,lambda)
  Heval->setValue(1, 2.0 - lambda(0));

#ifdef SICONOS_DEBUG
  std::cout << "modif heval : \n";
  Heval->display();
#endif

}



void adjointInput::computeg(double t, Interaction& inter)
{
  SiconosVector& lambda = *inter.lambda(0);
  SiconosVector workX = *inter.data(x);

#ifdef SICONOS_DEBUG
  std::cout << "************      computeG at: " << t << std::endl;
#endif

  SP::SiconosVector K2P(new SiconosVector(2));
  SP::SiconosVector P(new SiconosVector(2));
  P->setValue(0, workX(2));
  P->setValue(1, workX(3));



  prod(*K2, *P, *K2P, true);

  SP::SiconosVector betatmp(new SiconosVector(2));
  beta(t, workX, betatmp);

  (*inter.data(g_alpha)).setValue(0, betatmp->getValue(0) * (lambda(0) - 1.0));       //R=g_barre(x,lambda_barre)
  (*inter.data(g_alpha)).setValue(1, (betatmp->getValue(1)) * (lambda(0) - 1.0));
  (*inter.data(g_alpha)).setValue(2, (K2P->getValue(0)) * (lambda(0) - 1.0));
  (*inter.data(g_alpha)).setValue(3, (K2P->getValue(1)) * (lambda(0) - 1.0));


#ifdef SICONOS_DEBUG
  std::cout << "modif g_alpha : \n";
  data[g_alpha]->display();
#endif
}

/** default function to compute jacobianH
 *  \param double : current time
 *  \param index for jacobian (0: jacobian according to x, 1 according to lambda)
 */
void adjointInput::computeJachx(double t, Interaction& inter)
{
  SiconosVector workX = *inter.data(x);

  double *h = &(*_jachx)(0, 0);
#ifdef SICONOS_DEBUG
  std::cout << "computeJachx " << " at " << " " << t << std::endl;
#endif

  SP::SiconosVector betatmp(new SiconosVector(2));
  beta(t, workX, betatmp);

  SP::SiconosMatrix jacbetaXtmp(new SimpleMatrix(2, 2));

  JacobianXbeta(t, workX, jacbetaXtmp);


  h[0] = workX(3);
  h[2] = -workX(2) ;
  h[4] = (-workX(1) + 1.0) ;
  h[6] = workX(0);
  h[1] = 0.0;
  h[3] = 0.0 ;
  h[5] = 0.0 ;
  h[7] = 0.0;



#ifdef SICONOS_DEBUG
  std::cout << "modif Jachx : \n";
  _jachx->display();
#endif


}
void adjointInput::computeJachlambda(double t, Interaction& inter)
{
  double *h = &(*_jachlambda)(0, 0);
#ifdef SICONOS_DEBUG
  std::cout << "computeJachlambda " << " at " << " " << t << std::endl;
#endif

  h[0] = 0.0              ;
  h[2] = 1.0 ;
  h[1] = -1.0             ;
  h[3] = 0.0 ;

#ifdef SICONOS_DEBUG
  std::cout << "modif Jachlambda : \n";
  _jachlambda->display();
#endif

}
/** default function to compute jacobianG according to lambda
 *  \param double : current time
 *  \param index for jacobian: at the time only one possible jacobian => i = 0 is the default value .
 */
void adjointInput::computeJacgx(double t, Interaction& inter)
{

  double *g = &(*_jacgx)(0, 0);
#ifdef SICONOS_DEBUG
  std::cout << "computeJacgx " << " at " << " " << t << std::endl;
#endif

  SiconosVector& lambda = *inter.lambda(0);
  SiconosVector workX = *inter.data(x);

  SP::SiconosMatrix jacbetaXtmp(new SimpleMatrix(2, 2));

  JacobianXbeta(t, workX, jacbetaXtmp);

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
  _jacgx->display();
#endif


}
void adjointInput::computeJacglambda(double t, Interaction& inter)
{

  double *g = &(*_jacglambda)(0, 0);
  SiconosVector workX = *inter.data(x);
#ifdef SICONOS_DEBUG
  std::cout << "computeJacglambda " << " at " << " " << t << std::endl;
#endif

  SP::SiconosVector K2P(new SiconosVector(2));
  SP::SiconosVector P(new SiconosVector(2));
  P->setValue(0, workX(2));
  P->setValue(1, workX(3));

  prod(*K2, *P, *K2P, true);

  SP::SiconosVector betatmp(new SiconosVector(2));
  beta(t, workX, betatmp);

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
  _jacglambda->display();
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
