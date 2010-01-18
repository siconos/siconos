#ifndef ELECRELATION_CPP
#define ELECRELATION_CPP

#include "elecRelation.h"

#include "circuit.h"

elecRelation::elecRelation():
  FirstOrderType2R()
{
}


void elecRelation::initialize(SP::Interaction inter)
{
  FirstOrderType2R::initialize(inter);
  unsigned int sizeY = interaction()->getSizeOfY();
  unsigned int sizeDS = interaction()->getSizeOfDS();
  SP::SiconosVector y = interaction()->y(0);
  SP::SiconosVector lambda = interaction()->lambda(0);

  double t0 = 0;

  _workL.reset(new SimpleVector(interaction()->getSizeOfY()));
  Jachx->resize(sizeY, sizeDS);
  _jachlambda->resize(sizeY, sizeY);

  jacgx->resize(sizeDS, sizeDS);
  Jacglambda->resize(sizeDS, sizeY);

#ifdef CLSC_CIRCUIT
  _workX->setValue(0, 0);
  _workL->setValue(0, 20.0 / (sR2 + sR1s));
  _workL->setValue(1, -20.0 / (sR2 + sR1d));
  _workL->setValue(2, 20 - 20 * ((sR1s) / (sR2 + sR1s)));
  _workL->setValue(3, 0);
  _workL->setValue(4, sE_plus);
  _workL->setValue(5, 0);
  _workL->setValue(6, 0);
  _workL->setValue(7, 20 - 20 * ((sR1s) / (sR2 + sR1s)));
  _workL->setValue(8, sR2 - sR1d);
  y->setValue(0, 0);
  y->setValue(1, 0);
  y->setValue(2, 0);
  y->setValue(3, 0);
  y->setValue(4, 0);
  y->setValue(5, sR2 - sR1s);
  y->setValue(6, sE_plus);
  y->setValue(7, 0);
  y->setValue(8, 0);

#else
  _workX->setValue(0, 0);
  _workL->setValue(0, 0);
  _workL->setValue(1, 0);
  _workL->setValue(2, 0);
  _workL->setValue(3, 0);
#endif

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
void elecRelation::computeh(double t)
{

  *_workX = *data[x];
  SP::SiconosVector lambda = interaction()->lambda(0);
  *_workL = *lambda;

#ifdef SICONOS_DEBUGc
  std::cout << "********         computeh at " << t << std::endl;
#endif
  SP::SiconosVector Heval = interaction()->relation()->Halpha();
#ifdef CLSC_CIRCUIT
  Heval->setValue(0, _workL->getValue(4) - source(t));
  Heval->setValue(1, _workX->getValue(0) - (_workL->getValue(3)) / sR);
  Heval->setValue(2, _workL->getValue(2) - 20 + _workL->getValue(0) * (_workL->getValue(6) + sR1s));
  Heval->setValue(3, _workL->getValue(2) + _workL->getValue(1) * (_workL->getValue(8) + sR1d));
  Heval->setValue(4, _workX->getValue(0) - _workL->getValue(0) - _workL->getValue(1));
  Heval->setValue(5, sR2 - _workL->getValue(6) - sR1s);
  Heval->setValue(6, sAmpli * (_workL->getValue(4) - _workL->getValue(3)) + _workL->getValue(5));
  Heval->setValue(7, sR2 - _workL->getValue(8) - sR1d);
  Heval->setValue(8, -_workL->getValue(2) + _workL->getValue(7));
#else
  Heval->setValue(0, - _workL->getValue(0) + source(t));
  Heval->setValue(1, -_workL->getValue(0) + _workX->getValue(0) + (_workL->getValue(3) + sR1)* _workL->getValue(1));
  Heval->setValue(2, sR2 - _workL->getValue(3) - sR1);
  Heval->setValue(3, _workL->getValue(0) + _workL->getValue(2));
#endif

#ifdef SICONOS_DEBUG
  std::cout << "modif heval : \n";
  Heval->display();
#endif

}



void elecRelation::computeg(double t)
{
  SP::SiconosVector lambda = interaction()->lambda(0);
  *_workL = *lambda;

#ifdef SICONOS_DEBUG
  std::cout << "************      computeg at: " << t << std::endl;
#endif

#ifdef CLSC_CIRCUIT
  (*data[g_alpha]).setValue(0, (_workL->getValue(2) - _workL->getValue(3)) / sL);
#else
  (*data[g_alpha]).setValue(0, (_workL->getValue(1)) / sC);
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
void elecRelation::computeJachx(double t)
{

  SP::SiconosVector lambda = interaction()->lambda(0);
  *_workL = *lambda;
  double *h = &(*Jachx)(0, 0);
#ifdef SICONOS_DEBUG
  std::cout << "computeJachx " << " at " << " " << t << std::endl;
#endif

#ifdef CLSC_CIRCUIT
  h[0] = 0;
  h[1] = 1;
  h[2] = 0;
  h[3] = 0;
  h[4] = 1;
  h[5] = 0;
  h[6] = 0;
  h[7] = 0;
  h[8] = 0;
#else
  h[0] = 0;
  h[1] = 1;
  h[2] = 0;
  h[3] = 0;
#endif

}
void elecRelation::computeJachlambda(double t)
{

  SP::SiconosVector lambda = interaction()->lambda(0);
  *_workL = *lambda;
  double *h = &(*_jachlambda)(0, 0);
#ifdef SICONOS_DEBUG
  std::cout << "computeJachlambda " << " at " << " " << t << std::endl;
#endif

#ifdef CLSC_CIRCUIT
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
  h[2] = _workL->getValue(6) + sR1s;
  h[11] = 0;
  h[20] = 1;
  h[29] = 0;
  h[38] = 0;
  h[47] = 0;
  h[56] = _workL->getValue(0);
  h[65] = 0;
  h[74] = 0;
  h[3] = 0;
  h[12] = _workL->getValue(8) + sR1d;
  h[21] = 1;
  h[30] = 0;
  h[39] = 0;
  h[48] = 0;
  h[57] = 0;
  h[66] = 0;
  h[75] = _workL->getValue(1);
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
#else
  h[0] = -1;
  h[4] = 0;
  h[8] = 0;
  h[12] = 0;
  h[1] = -1;
  h[5] = _workL->getValue(3) + sR1;
  h[9] = 0;
  h[13] = _workL->getValue(1);
  h[2] = 0;
  h[6] = 0;
  h[10] = 0;
  h[14] = -1;
  h[3] = 1;
  h[7] = 0;
  h[11] = 1;
  h[15] = 0;
#endif

}
/** default function to compute jacobianG according to lambda
 *  \param double : current time
 *  \param index for jacobian: at the time only one possible jacobian => i = 0 is the default value .
 */
void elecRelation::computeJacgx(double t)
{

  double *g = &(*jacgx)(0, 0);
#ifdef SICONOS_DEBUG
  std::cout << "computeJacgx " << " at " << " " << t << std::endl;
#endif
#ifdef CLSC_CIRCUIT
  g[0] = 0;
#else
  g[0] = 0;
#endif
}
void elecRelation::computeJacglambda(double t)
{

  double *g = &(*Jacglambda)(0, 0);
#ifdef SICONOS_DEBUG
  std::cout << "computeJacglambda " << " at " << " " << t << std::endl;
#endif
#ifdef CLSC_CIRCUIT
  g[0] = 0;
  g[1] = 0;
  g[2] = 1 / sL;
  g[3] = -1 / sL;
  g[4] = 0;
  g[5] = 0;
  g[6] = 0;
  g[7] = 0;
  g[8] = 0;
#else
  g[0] = 0;
  g[1] = 1 / sC;
  g[2] = 0;
  g[3] = 0;
#endif
}
#endif
