#include "mlcp_enum_tool.h"
#include <stdio.h>
#include "NumericsOptions.h"

static unsigned long long int sCurrentEnum = 0;
static unsigned long long int sCmpEnum = 0;
static unsigned long long int sNbCase = 0;
static double sProgress = 0;
static int sMm = 0;


void initEnum(int M)
{
  int cmp;
  /*  sCurrentEnum = 0;*/

  sCmpEnum = 0;
  sNbCase = 1;
  sMm = M;

  for (cmp = 0; cmp < sMm; cmp++)
    sNbCase = sNbCase << 1;
  sProgress = 0;
}

void affectW2V(int * W2V)
{
  unsigned long  int aux = sCurrentEnum;
  for (unsigned int i = 0; i < sMm; i++)
  {
    W2V[i] = aux & 1;
    aux = aux >> 1;
  }
  if (verbose)
  {
    for (unsigned int i = 0; i < sMm; i++)
      printf("wv[%d]=%d \t", i, W2V[i]);
    printf("\n");
  }

}

int nextEnum(int * W2V)
{
  if (sCmpEnum == sNbCase)
    return 0;
  if (sCurrentEnum >= sNbCase)
  {
    sCurrentEnum = 0;
  }
  if (verbose)
    printf("try enum :%d\n", (int)sCurrentEnum);
  affectW2V(W2V);
  sCurrentEnum++;
  sCmpEnum++;
  if (verbose && sCmpEnum > sProgress * sNbCase)
  {
    sProgress += 0.001;
    printf(" progress %f %d \n", sProgress, (int) sCurrentEnum);
  }

  return 1;
}
