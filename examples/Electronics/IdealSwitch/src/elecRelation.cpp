#ifndef ELECRELATION_CPP
#define ELECRELATION_CPP

#include "elecRelation.h"

#include "circuit.h"

elecRelation::elecRelation():
  FirstOrderType2R()
{
}


void elecRelation::initializeWorkVectorsAndMatrices(Interaction& inter, VectorOfBlockVectors& DSlink, VectorOfVectors& workV, VectorOfSMatrices& workM)
{
  FirstOrderType2R::initializeWorkVectorsAndMatrices(inter, DSlink, workV, workM);
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
void elecRelation::computeh(double t, SiconosVector& x, SiconosVector& lambda, SiconosVector& y)
{


#ifdef CLSC_CIRCUIT
  y(0) = lambda(4) - source(t);
  y(1) = x(0) - (lambda(3)) / sR;
  y(2) = lambda(2) - 20 + lambda(0) * (lambda(6) + sR1s);
  y(3) = lambda(2) + lambda(1) * (lambda(8) + sR1d);
  y(4) = x(0) - lambda(0) - lambda(1);
  y(5) = sR2 - lambda(6) - sR1s;
  y(6) = sAmpli * (lambda(4) - lambda(3)) + lambda(5);
  y(7) = sR2 - lambda(8) - sR1d;
  y(8) = -lambda(2) + lambda(7);
#else
  y(0) = - lambda(0) + source(t);
  y(1) = -lambda(0) + x(0) + (lambda(3) + sR1)* lambda(1);
  y(2) = sR2 - lambda(3) - sR1;
  y(3) = lambda(0) + lambda(2);
#endif

}



void elecRelation::computeg(double t, SiconosVector& lambda, SiconosVector& r)
{
#ifdef SICONOS_DEBUG
  std::cout << "************      computeg at: " << t << std::endl;
#endif

#ifdef CLSC_CIRCUIT
  r(0) = (lambda(2) - lambda(3)) / sL;
#else
  r(0) = lambda(1) / sC;
#endif

#ifdef SICONOS_DEBUG
  std::cout << "modif g_alpha : \n";
  r.display();
#endif
}

/** default function to compute jacobianH
 *  \param double : current time
 *  \param index for jacobian (0: jacobian according to x, 1 according to lambda)
 */
void elecRelation::computeJachx(double t, SiconosVector& x, SiconosVector& lambda, SimpleMatrix& C)
{

  double* h = C.getArray();
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
void elecRelation::computeJachlambda(double t, SiconosVector& x, SiconosVector& lambda, SimpleMatrix& D)
{

  double* h = D.getArray();
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
  h[2] = lambda(6) + sR1s;
  h[11] = 0;
  h[20] = 1;
  h[29] = 0;
  h[38] = 0;
  h[47] = 0;
  h[56] = lambda(0);
  h[65] = 0;
  h[74] = 0;
  h[3] = 0;
  h[12] = lambda(8) + sR1d;
  h[21] = 1;
  h[30] = 0;
  h[39] = 0;
  h[48] = 0;
  h[57] = 0;
  h[66] = 0;
  h[75] = lambda(1);
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
  h[5] = lambda(3) + sR1;
  h[9] = 0;
  h[13] = lambda(1);
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

void elecRelation::computeJacglambda(double time, SiconosVector& lambda, SimpleMatrix& B)
{

  double *g = B.getArray();
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
