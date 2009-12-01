#ifndef ADJOINTINPUT_CPP
#define ADJOINTINPUT_CPP

#include "adjointInput.h"

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

  double t0 = 0;

  workL.reset(new SimpleVector(interaction()->getSizeOfY()));
  JacXH->resize(sizeY, sizeDS);
  JacLH->resize(sizeY, sizeY);

  JacXG->resize(sizeDS, sizeDS);
  JacLG->resize(sizeDS, sizeY);


  workX->setValue(0, 0);
  workL->setValue(0, 0);
  workL->setValue(1, 0);

  *lambda = *workL;

  //  computeH(t0);
  computeG(t0);
  computeJacH(t0);
  computeJacG(t0);
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
void adjointInput::computeH(double t)
{

  *workX = *data[x];
  SP::SiconosVector lambda = interaction()->lambda(0);
  *workL = *lambda;

#ifdef SICONOS_DEBUGc
  std::cout << "********         computeH at " << t << std::endl;
#endif
  SP::SiconosVector Heval = interaction()->relation()->Halpha();

  double betap = (1.0 / 2.0 * workX->getValue(1) + 1.0 / 2.0) * workX->getValue(2) +
                 (-1.0 / 2.0 * workX->getValue(0)) * workX->getValue(3);


  Heval->setValue(0, betap + workL->getValue(1));
  Heval->setValue(1, 2.0 - workL->getValue(0));

#ifdef SICONOS_DEBUG
  std::cout << "modif heval : \n";
  Heval->display();
#endif

}



void adjointInput::computeG(double t)
{
  SP::SiconosVector lambda = interaction()->lambda(0);
  *workL = *lambda;
  *workX = *data[x];

#ifdef SICONOS_DEBUG
  std::cout << "************      computeG at: " << t << std::endl;
#endif



  SP::SimpleMatrix K2(new SimpleMatrix(2, 2));
  K2->setValue(0, 0, 0.0);
  K2->setValue(0, 1, +1.0 / 2.0);
  K2->setValue(1, 0, -1.0 / 2.0);
  K2->setValue(1, 1, 0.0);

  SP::SimpleVector K2P(new SimpleVector(2));
  SP::SimpleVector P(new SimpleVector(2));
  P->setValue(0, workX->getValue(2));
  P->setValue(1, workX->getValue(3));



  prod(*K2, *P, *K2P, true);

  (*data[g_alpha]).setValue(0, -(1.0 / 2.0 * workX->getValue(1) + 1.0 / 2.0) * (workL->getValue(0)));
  (*data[g_alpha]).setValue(1, -(-1.0 / 2.0 * workX->getValue(0)) * (workL->getValue(0)));
  (*data[g_alpha]).setValue(2, -(K2P->getValue(0)) * (workL->getValue(1)));
  (*data[g_alpha]).setValue(3, -(K2P->getValue(1)) * (workL->getValue(1)));


#ifdef SICONOS_DEBUG
  std::cout << "modif g_alpha : \n";
  data[g_alpha]->display();
#endif
}

/** default function to compute jacobianH
 *  \param double : current time
 *  \param index for jacobian (0: jacobian according to x, 1 according to lambda)
 */
void adjointInput::computeJacXH(double t)
{

  SP::SiconosVector lambda = interaction()->lambda(0);
  *workL = *lambda;
  *workX = *data[x];

  double *h = &(*JacXH)(0, 0);
#ifdef SICONOS_DEBUG
  std::cout << "computeJacXH " << " at " << " " << t << std::endl;
#endif


  h[0] = -1.0 / 2.0 * workX->getValue(3);
  h[2] = 1.0 / 2.0 * workX->getValue(2);
  h[4] = 1.0 / 2.0 * (workX->getValue(1) + 1.0);
  h[6] = -1.0 / 2.0 * workX->getValue(0);
  h[1] = 0.0;
  h[3] = 0.0;
  h[5] = 0.0;
  h[7] = 0.0;


}
void adjointInput::computeJacLH(double t)
{

  SP::SiconosVector lambda = interaction()->lambda(0);
  *workL = *lambda;
  double *h = &(*JacLH)(0, 0);
#ifdef SICONOS_DEBUG
  std::cout << "computeJacLH " << " at " << " " << t << std::endl;
#endif


  h[0] = 0.0;
  h[2] = -1.0;
  h[2] = -1.0;
  h[3] = 0.0;


}
/** default function to compute jacobianG according to lambda
 *  \param double : current time
 *  \param index for jacobian: at the time only one possible jacobian => i = 0 is the default value .
 */
void adjointInput::computeJacXG(double t)
{

  double *g = &(*JacXG)(0, 0);
#ifdef SICONOS_DEBUG
  std::cout << "computeJacXG " << " at " << " " << t << std::endl;
#endif

  SP::SiconosVector lambda = interaction()->lambda(0);
  *workL = *lambda;
  *workX = *data[x];


  SP::SimpleMatrix K2(new SimpleMatrix(2, 2));
  K2->setValue(0, 0, 0.0);
  K2->setValue(0, 1, +1.0 / 2.0);
  K2->setValue(1, 0, -1.0 / 2.0);
  K2->setValue(1, 1, 0.0);


  g[0] = 0.0;
  g[4] = -1.0 / 2.0 * workL->getValue(0)  ;
  g[8] = 0.0;
  g[12] = 0.0;
  g[1] = +1.0 / 2.0 * workL->getValue(1) ;
  g[5] = 0.0 ;
  g[9] = 0.0  ;
  g[13] = 0.0 ;
  g[2] = 0.0;
  g[6] = 0.0;
  g[10] = -K2->getValue(0, 0) * workL->getValue(0);
  g[14] = -K2->getValue(0, 1) * workL->getValue(0);
  g[3] = 0.0;
  g[7] = 0.0;
  g[11] = -K2->getValue(1, 0) * workL->getValue(1);
  g[15] = -K2->getValue(1, 1) * workL->getValue(1);




}
void adjointInput::computeJacLG(double t)
{

  double *g = &(*JacLG)(0, 0);
#ifdef SICONOS_DEBUG
  std::cout << "computeJacLG " << " at " << " " << t << std::endl;
#endif
  SP::SimpleMatrix K2(new SimpleMatrix(2, 2));
  K2->setValue(0, 0, 0.0);
  K2->setValue(0, 1, +1.0 / 2.0);
  K2->setValue(1, 0, -1.0 / 2.0);
  K2->setValue(1, 1, 0.0);
  *workX = *data[x];

  SP::SimpleVector K2P(new SimpleVector(2));
  SP::SimpleVector P(new SimpleVector(2));
  P->setValue(0, workX->getValue(2));
  P->setValue(1, workX->getValue(3));



  prod(*K2, *P, *K2P, true);


  g[0] = -1.0 / 2.0 * (workX->getValue(1) + 1.0);
  g[4] = 0.0;
  g[1] = 0.0 ;
  g[5] = 1.0 / 2.0 * (workX->getValue(0));
  g[2] = -K2P->getValue(0);
  g[6] = 0.0;
  g[3] = 0.0;
  g[7] = -K2P->getValue(1) ;



}
#endif
