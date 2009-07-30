#ifndef ELECRELATION_CPP
#define ELECRELATION_CPP

#include "elecRelation.h"

#include "circuit.h"

elecRelation::elecRelation():
  FirstOrderType2R()
{
  JacH.resize(2);
  JacG.resize(2);
  JacH[0].reset(new PluggedMatrix());
  JacH[1].reset(new PluggedMatrix());
  JacG[0].reset(new PluggedMatrix());
  JacG[1].reset(new PluggedMatrix());
}


void elecRelation::initialize(SP::Interaction inter)
{
  FirstOrderType2R::initialize(inter);
  unsigned int sizeY = getInteractionPtr()->getSizeOfY();
  unsigned int sizeDS = getInteractionPtr()->getSizeOfDS();
  SP::SiconosVector y = getInteractionPtr()->getYPtr(0);
  SP::SiconosVector lambda = getInteractionPtr()->getLambdaPtr(0);

  double t0 = 0;

  workL.reset(new SimpleVector(getInteractionPtr()->getSizeOfY()));
  JacH[0]->resize(sizeY, sizeDS);
  JacH[1]->resize(sizeY, sizeY);

  JacG[0]->resize(sizeDS, sizeDS);
  JacG[1]->resize(sizeDS, sizeY);

#ifdef CLSC_CIRCUIT
  workX->setValue(0, 0);
  workL->setValue(0, 20.0 / (sR2 + sR1));
  workL->setValue(1, -20.0 / (sR2 + sR1));
  workL->setValue(2, 20 - 20 * ((sR1) / (sR2 + sR1)));
  workL->setValue(3, 0);
  workL->setValue(4, sE_plus);
  workL->setValue(5, 0);
  workL->setValue(6, 0);
  workL->setValue(7, 20 - 20 * ((sR1) / (sR2 + sR1)));
  workL->setValue(8, sR2 - sR1);
  y->setValue(0, 0);
  y->setValue(1, 0);
  y->setValue(2, 0);
  y->setValue(3, 0);
  y->setValue(4, 0);
  y->setValue(5, sR2 - sR1);
  y->setValue(6, sE_plus);
  y->setValue(7, 0);
  y->setValue(8, 0);

#else
  workX->setValue(0, 0);
  workL->setValue(0, 0);
  workL->setValue(1, 0);
  workL->setValue(2, 0);
  workL->setValue(3, 0);
#endif

  *lambda = *workL;

  //  computeH(t0);
  computeG(t0);
  computeJacH(t0, 0);
  computeJacH(t0, 1);
  computeJacG(t0, 0);
  computeJacG(t0, 1);
  *data[r] = *data[g_alpha];
#ifdef SICONOS_DEBUG
  std::cout << "data[r (g_alpha)] init\n";
  data[r]->display();
#endif

}


double elecRelation::source(double t)
{
  double daux = 0;
#ifdef CLSC_CIRCUIT
  double numT = t / sT;
  int aux = (int) floor(numT);
  daux = sE_plus - ((sE_plus - sE_moins) / sT) * t + (sE_plus - sE_moins) * aux;
#ifdef SICONOS_DEBUG
  std::cout << "source(" << t << ")=" << daux << std::endl;
#endif
  return daux;

#else
  daux =  sin(sW * t);
#ifdef SICONOS_DEBUG
  std::cout << "source(" << t << ")=" << daux << std::endl;
#endif
  return daux;
#endif
}

/*y = h(X,lambda)*/
void elecRelation::computeH(double t)
{

  *workX = *data[x];
  SP::SiconosVector lambda = getInteractionPtr()->getLambdaPtr(0);
  *workL = *lambda;

#ifdef SICONOS_DEBUG
  std::cout << "********         computeH at " << t << std::endl;
#endif
  SP::SiconosVector Heval = getInteractionPtr()->getRelationPtr()->getHalphaPtr();
#ifdef CLSC_CIRCUIT
  Heval->setValue(0, workL->getValue(4) - source(t));
  Heval->setValue(1, workX->getValue(0) - (workL->getValue(3)) / sR);
  Heval->setValue(2, workL->getValue(2) - 20 + workL->getValue(0) * (workL->getValue(6) + sR1));
  Heval->setValue(3, workL->getValue(2) + workL->getValue(1) * (workL->getValue(8) + sR1));
  Heval->setValue(4, workX->getValue(0) - workL->getValue(0) - workL->getValue(1));
  Heval->setValue(5, sR2 - workL->getValue(6) - sR1);
  Heval->setValue(6, sAmpli * (workL->getValue(4) - workL->getValue(3)) + workL->getValue(5));
  Heval->setValue(7, sR2 - workL->getValue(8) - sR1);
  Heval->setValue(8, -workL->getValue(2) + workL->getValue(7));
#else
  Heval->setValue(0, - workL->getValue(0) + source(t));
  Heval->setValue(1, -workL->getValue(0) + workX->getValue(0) + (workL->getValue(3) + sR1)* workL->getValue(1));
  Heval->setValue(2, sR2 - workL->getValue(3) - sR1);
  Heval->setValue(3, workL->getValue(0) + workL->getValue(2));
#endif

#ifdef SICONOS_DEBUG
  std::cout << "modif heval : \n";
  Heval->display();
#endif

}



void elecRelation::computeG(double t)
{
  SP::SiconosVector lambda = getInteractionPtr()->getLambdaPtr(0);
  *workL = *lambda;

#ifdef SICONOS_DEBUG
  std::cout << "************      computeG at: " << t << std::endl;
#endif

#ifdef CLSC_CIRCUIT
  (*data[g_alpha]).setValue(0, (workL->getValue(2) - workL->getValue(3)) / sL);
#else
  (*data[g_alpha]).setValue(0, (workL->getValue(1)) / sC);
#endif

#ifdef SICONOS_DEBUG
  std::cout << "modif g_alpha : \n";
  data[g_alpha]->display();
#endif
}

/** default function to compute jacobianH
 *  \param double : current time
 *  \param index for jacobian (0: jacobian according to x, 1 according to lambda)
 */
void elecRelation::computeJacH(double t, unsigned int index)
{

  SP::SiconosVector lambda = getInteractionPtr()->getLambdaPtr(0);
  *workL = *lambda;
  double *h = &(*(JacH[index]))(0, 0);
#ifdef SICONOS_DEBUG
  std::cout << "computeJacH " << index << " at " << " " << t << std::endl;
#endif

#ifdef CLSC_CIRCUIT
  if (index == 0)
  {
    h[0] = 0;
    h[1] = 1;
    h[2] = 0;
    h[3] = 0;
    h[4] = 1;
    h[5] = 0;
    h[6] = 0;
    h[7] = 0;
    h[8] = 0;
  }
  else
  {
    h[0] = 0;
    h[9] = 0;
    h[18] = 0;
    h[27] = 0;
    h[36] = 1;
    h[45] = 0;
    h[54] = 0;
    h[63] = 0;
    h[72] = 0;
    h[1] = 0;
    h[10] = 0;
    h[19] = 0;
    h[28] = -1 / sR;
    h[37] = 0;
    h[46] = 0;
    h[55] = 0;
    h[64] = 0;
    h[73] = 0;
    h[2] = workL->getValue(6) + sR1;
    h[11] = 0;
    h[20] = 1;
    h[29] = 0;
    h[38] = 0;
    h[47] = 0;
    h[56] = workL->getValue(0);
    h[65] = 0;
    h[74] = 0;
    h[3] = 0;
    h[12] = workL->getValue(8) + sR1;
    h[21] = 1;
    h[30] = 0;
    h[39] = 0;
    h[48] = 0;
    h[57] = 0;
    h[66] = 0;
    h[75] = workL->getValue(1);
    h[4] = -1;
    h[13] = -1;
    h[22] = 0;
    h[31] = 0;
    h[40] = 0;
    h[49] = 0;
    h[58] = 0;
    h[67] = 0;
    h[76] = 0;
    h[5] = 0;
    h[14] = 0;
    h[23] = 0;
    h[32] = 0;
    h[41] = 0;
    h[50] = 0;
    h[59] = -1;
    h[68] = 0;
    h[77] = 0;
    h[6] = 0;
    h[15] = 0;
    h[24] = 0;
    h[33] = -sAmpli;
    h[42] = sAmpli;
    h[51] = 1;
    h[60] = 0;
    h[69] = 0;
    h[78] = 0;
    h[7] = 0;
    h[16] = 0;
    h[25] = 0;
    h[34] = 0;
    h[43] = 0;
    h[52] = 0;
    h[61] = 0;
    h[70] = 0;
    h[79] = -1;
    h[8] = 0;
    h[17] = 0;
    h[26] = -1;
    h[35] = 0;
    h[44] = 0;
    h[53] = 0;
    h[62] = 0;
    h[71] = 1;
    h[80] = 0;
  }
#else
  if (index == 0)
  {
    h[0] = 0;
    h[1] = 1;
    h[2] = 0;
    h[3] = 0;
  }
  else
  {
    h[0] = -1;
    h[4] = 0;
    h[8] = 0;
    h[12] = 0;
    h[1] = -1;
    h[5] = workL->getValue(3) + sR1;
    h[9] = 0;
    h[13] = workL->getValue(1);
    h[2] = 0;
    h[6] = 0;
    h[10] = 0;
    h[14] = -1;
    h[3] = 1;
    h[7] = 0;
    h[11] = 1;
    h[15] = 0;
  }
#endif

}

/** default function to compute jacobianG according to lambda
 *  \param double : current time
 *  \param index for jacobian: at the time only one possible jacobian => i = 0 is the default value .
 */
void elecRelation::computeJacG(double t, unsigned int index)
{

  double *g = &(*(JacG[index]))(0, 0);
#ifdef SICONOS_DEBUG
  std::cout << "computeJacG " << index << " at " << " " << t << std::endl;
#endif
#ifdef CLSC_CIRCUIT
  if (index == 0)
  {
    g[0] = 0;
  }
  else
  {
    g[0] = 0;
    g[1] = 0;
    g[2] = 1 / sL;
    g[3] = -1 / sL;
    g[4] = 0;
    g[5] = 0;
    g[6] = 0;
    g[7] = 0;
    g[8] = 0;
  }
#else
  if (index == 0)
  {
    g[0] = 0;
  }
  else
  {
    g[0] = 0;
    g[1] = 1 / sC;
    g[2] = 0;
    g[3] = 0;
  }
#endif
}

#endif
