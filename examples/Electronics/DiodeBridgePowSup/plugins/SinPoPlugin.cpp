#ifdef _WIN32 
#define SICONOS_EXPORT extern "C" __declspec(dllexport) 
#else 
#define SICONOS_EXPORT extern "C" 
#endif  
#include <stdio.h>
#include <math.h>

SICONOS_EXPORT void SinPo(double time, unsigned int sizeOfE, double* EPtr, unsigned int sizeOfZ, double* z)
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
