#include <stdio.h>
#include <math.h>

extern "C"   void SinPo(double time, unsigned int sizeOfE, double* EPtr, unsigned int sizeOfZ, double* z)
{

  double omega = 1e4;
  double Voffset = 0.0;
  double amplitude = 10.0;
  double phase = 0.0;
  double VSinPo;

  VSinPo = Voffset + (amplitude * cos((omega * time) + phase));

  for (unsigned int i = 0; i < sizeOfE; i++) EPtr[i] = z[i];
  EPtr[2] -= VSinPo;
  EPtr[3] += VSinPo;
  z[4] = VSinPo;
}
