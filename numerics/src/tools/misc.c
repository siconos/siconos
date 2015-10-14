#include <stdio.h>
#include "misc.h"

void printm(unsigned int nl, unsigned int nc, double *m)
{
  for (unsigned int i = 0; i < nl; ++i)
  {
    for (unsigned int j = 0; j < nc; ++j)
    {
      printf("%32.24e ", *(m + j + nc * i));
    }
    printf("\n");
  }
}
