#include "enum_tool.h"
#include "numerics_verbose.h"                   // for verbose
#include <stdio.h>                              // for printf
#include <stdlib.h>                              // for printf
unsigned long long int enum_compute_nb_cases(int M)
{
  unsigned long long int nbCase = 1;
  for(int cmp = 0; cmp < M; cmp++)
    nbCase = nbCase << 1;
  return nbCase;
}

EnumerationStruct * enum_init(int M)
{
  EnumerationStruct * enum_struct = (EnumerationStruct *)malloc(sizeof(EnumerationStruct));

  enum_struct->current = 0;
  /*  scurrent = 0;*/
  numerics_printf_verbose(1,"----- enum_init -- problem size :%i", M);
  numerics_printf_verbose(1,"----- enum_init -- currentEnum :%i", (int) enum_struct->current);

  enum_struct->counter = 0;

  enum_struct->nb_cases  = enum_compute_nb_cases(M);

  enum_struct->progress = 0;
  numerics_printf_verbose(1,"----- initEnum -- number of cases :%i", (int)enum_struct->nb_cases);
  return enum_struct;
}

static void enum_affect_zw(int * zw, int size, EnumerationStruct * enum_struct)
{
  unsigned long  int aux = enum_struct->current;
  for(int i = 0; i < size; i++)
  {
    zw[i] = aux & 1;
    aux = aux >> 1;
  }

  if(verbose > 1)
  {
    for(int i = 0; i < size; i++)
      printf("zw[%d]=%d \t", i, zw[i]);
    printf("\n");
  }
}
int enum_next(int * zw, int size, EnumerationStruct * enum_struct)
{
  if(enum_struct->counter == enum_struct->nb_cases)
    return 0;
  if(enum_struct->current >= enum_struct->nb_cases)
    enum_struct->current = 0;

  numerics_printf_verbose(1,"----- enum_next -- try enum :%d", (int)enum_struct->current);

  enum_affect_zw(zw, size, enum_struct);
  enum_struct->current++;
  enum_struct->counter++;

  if(verbose && enum_struct->counter > (unsigned long int) enum_struct->progress * enum_struct->nb_cases)
  {
    enum_struct->progress  += 0.001;
    numerics_printf_verbose(1,"progress %f %d / %d", enum_struct->progress, (int)enum_struct->current,  enum_struct->nb_cases);
  }

  return 1;
}
