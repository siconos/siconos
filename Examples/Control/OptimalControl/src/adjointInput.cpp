#ifndef ADJOINTINPUT_CPP
#define ADJOINTINPUT_CPP

#include "adjointInput.h"
//#define SICONOS_DEBUG
adjointInput::adjointInput():
  FirstOrderType2R()
{
}


void adjointInput::initialize(SP::Interaction inter)
{
  FirstOrderType2R::initialize(inter);
  unsigned int sizeY = interaction()->getSizeOfY();
  unsigned int sizeDS = interaction()->getSizeOfDS();
  SP::SiconosVector y = interaction()->y(0);
  SP::SiconosVector lambda = interaction()->lambda(0);

  K2 = new SimpleMatrix(2, 2);
  K2->setValue(0, 0, 0.0);
  K2->setValue(0, 1, -1.0 / 2.0);
  K2->setValue(1, 0, 1.0 / 2.0);
  K2->setValue(1, 1, 0.0);

  double t0 = 0;

  _workL.reset(new SiconosVector(interaction()->getSizeOfY()));
  Jachx->resize(sizeY, sizeDS);
  _jachlambda->resize(sizeY, sizeY);

  jacgx->resize(sizeDS, sizeDS);
  Jacglambda->resize(sizeDS, sizeY);


  _workX->setValue(0, 0);
  _workL->setValue(0, 0);
  _workL->setValue(1, 0);

  *lambda = *_workL;

  //  computeH(t0);
  computeg(t0);
  computeJach(t0);
  computeJacg(t0);
  *data[r] = *data[g_alpha];
#ifdef SICONOS_DEBUG
  std::cout << "data[r (g_alpha)] init\n";
  data[r]->display();
#endif

}


double adjointInput::source(double t)
{
  double daux = 0;
  return daux;
}

/*y = h(X,lambda)*/
void adjointInput::computeh(double t)
{

  *_workX = *data[x];
  SP::SiconosVector lambda = interaction()->lambda(0);
  *_workL = *lambda;

#ifdef SICONOS_DEBUG
  std::cout << "********         computeH at " << t << std::endl;
#endif
  SP::SiconosVector Heval = interaction()->relation()->Halpha();


  SP::SiconosVector betatmp(new SiconosVector(2));
  beta(t, _workX, betatmp);

  double betap = 2.0 * ((betatmp->getValue(0)) * _workX->getValue(2) + (betatmp->getValue(1)) * _workX->getValue(3));



  Heval->setValue(0, _workL->getValue(1) + betap);   //y_barre =Heval(x,lambda)
  Heval->setValue(1, 2.0 - _workL->getValue(0));

#ifdef SICONOS_DEBUG
  std::cout << "modif heval : \n";
  Heval->display();
#endif

}



void adjointInput::computeg(double t)
{
  SP::SiconosVector lambda = interaction()->lambda(0);
  *_workL = *lambda;
  *_workX = *data[x];

#ifdef SICONOS_DEBUG
  std::cout << "************      computeG at: " << t << std::endl;
#endif

  SP::SiconosVector K2P(new SiconosVector(2));
  SP::SiconosVector P(new SiconosVector(2));
  P->setValue(0, _workX->getValue(2));
  P->setValue(1, _workX->getValue(3));



  prod(*K2, *P, *K2P, true);

  SP::SiconosVector betatmp(new SiconosVector(2));
  beta(t, _workX, betatmp);

  (*data[g_alpha]).setValue(0, betatmp->getValue(0) * (_workL->getValue(0) - 1.0));       //R=g_barre(x,lambda_barre)
  (*data[g_alpha]).setValue(1, (betatmp->getValue(1)) * (_workL->getValue(0) - 1.0));
  (*data[g_alpha]).setValue(2, (K2P->getValue(0)) * (_workL->getValue(0) - 1.0));
  (*data[g_alpha]).setValue(3, (K2P->getValue(1)) * (_workL->getValue(0) - 1.0));


#ifdef SICONOS_DEBUG
  std::cout << "modif g_alpha : \n";
  data[g_alpha]->display();
#endif
}

/** default function to compute jacobianH
 *  \param double : current time
 *  \param index for jacobian (0: jacobian according to x, 1 according to lambda)
 */
void adjointInput::computeJachx(double t)
{

  SP::SiconosVector lambda = interaction()->lambda(0);
  *_workL = *lambda;
  *_workX = *data[x];

  double *h = &(*Jachx)(0, 0);
#ifdef SICONOS_DEBUG
  std::cout << "computeJachx " << " at " << " " << t << std::endl;
#endif

  SP::SiconosVector betatmp(new SiconosVector(2));
  beta(t, _workX, betatmp);

  SP::SiconosMatrix jacbetaXtmp(new SimpleMatrix(2, 2));

  JacobianXbeta(t, _workX, jacbetaXtmp);


  h[0] = _workX->getValue(3);
  h[2] = -_workX->getValue(2) ;
  h[4] = (-_workX->getValue(1) + 1.0) ;
  h[6] = _workX->getValue(0);
  h[1] = 0.0;
  h[3] = 0.0 ;
  h[5] = 0.0 ;
  h[7] = 0.0;



#ifdef SICONOS_DEBUG
  std::cout << "modif Jachx : \n";
  Jachx->display();
#endif


}
void adjointInput::computeJachlambda(double t)
{

  SP::SiconosVector lambda = interaction()->lambda(0);
  *_workL = *lambda;
  *_workX = *data[x];
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
  Jachlambda->display();
#endif

}
/** default function to compute jacobianG according to lambda
 *  \param double : current time
 *  \param index for jacobian: at the time only one possible jacobian => i = 0 is the default value .
 */
void adjointInput::computeJacgx(double t)
{

  double *g = &(*jacgx)(0, 0);
#ifdef SICONOS_DEBUG
  std::cout << "computeJacgx " << " at " << " " << t << std::endl;
#endif

  SP::SiconosVector lambda = interaction()->lambda(0);
  *_workL = *lambda;
  *_workX = *data[x];

  SP::SiconosMatrix jacbetaXtmp(new SimpleMatrix(2, 2));

  JacobianXbeta(t, _workX, jacbetaXtmp);

  g[0] = jacbetaXtmp->getValue(0, 0) * (_workL->getValue(0) - 1.0) ;
  g[4] = jacbetaXtmp->getValue(0, 1) * (_workL->getValue(0) - 1.0)  ;
  g[8] = 0.0;
  g[12] = 0.0;
  g[1] = jacbetaXtmp->getValue(1, 0) * (_workL->getValue(0) - 1.0) ;
  g[5] = jacbetaXtmp->getValue(1, 1) * (_workL->getValue(0) - 1.0)  ;
  g[9] = 0.0  ;
  g[13] = 0.0 ;
  g[2] = 0.0;
  g[6] = 0.0;
  g[10] = K2->getValue(0, 0) * (_workL->getValue(0) - 1.0);
  g[14] = K2->getValue(0, 1) * (_workL->getValue(0) - 1.0);
  g[3] = 0.0;
  g[7] = 0.0;
  g[11] = K2->getValue(1, 0) * (_workL->getValue(0) - 1.0);
  g[15] = K2->getValue(1, 1) * (_workL->getValue(0) - 1.0);

#ifdef SICONOS_DEBUG
  std::cout << "modif Jacgx : \n";
  Jacgx->display();
#endif


}
void adjointInput::computeJacglambda(double t)
{

  double *g = &(*Jacglambda)(0, 0);
#ifdef SICONOS_DEBUG
  std::cout << "computeJacglambda " << " at " << " " << t << std::endl;
#endif

  SP::SiconosVector K2P(new SiconosVector(2));
  SP::SiconosVector P(new SiconosVector(2));
  P->setValue(0, _workX->getValue(2));
  P->setValue(1, _workX->getValue(3));

  prod(*K2, *P, *K2P, true);

  SP::SiconosVector betatmp(new SiconosVector(2));
  beta(t, _workX, betatmp);

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
  Jacglambda->display();
#endif

}


void adjointInput::beta(double t, SP::SiconosVector xvalue, SP::SiconosVector beta)
{


  beta->setValue(0, -1.0 / 2.0 * xvalue->getValue(1) + 1.0 / 2.0) ;
  beta->setValue(1, 1.0 / 2.0 * xvalue->getValue(0)) ;
#ifdef SICONOS_DEBUG
  std::cout << "beta\n" << std::endl;;
  beta->display();
#endif

}

void adjointInput::JacobianXbeta(double t, SP::SiconosVector xvalue, SP::SiconosMatrix JacXbeta)
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
