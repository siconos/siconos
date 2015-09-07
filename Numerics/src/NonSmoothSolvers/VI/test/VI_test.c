#include "VariationalInequality.h"
#include "stdlib.h"

#pragma GCC diagnostic ignored "-Wmissing-prototypes"

void Ftest(void * viIn, int n, double *x, double *F)
{
  int i;
  VariationalInequality * vi = (VariationalInequality* ) viIn;
  printf("Size of vi :%i\n", vi->size);
  
  for (i =0; i< vi->size ; i++)
  {
    F[i] = x[i];
  }
}
void PXtest(void *viIn, double *x, double *PX)
{
  VariationalInequality * vi = (VariationalInequality* ) viIn;
  printf("Size of vi :%i\n", vi->size);
  int i;
  for (i =0; i< vi->size ; i++)
  {
    PX[i] = x[i];
    if (PX[i] <0) PX[i]=0.0;
  }
}




int main(void)
{

  VariationalInequality vi;
  variationalInequality_clear(&vi);
  
  vi.size=10;
  //vi.Callback = (CallbackVI *)malloc(sizeof(CallbackVI));
  
  vi.env = &vi;
  
  vi.F = &Ftest;
  vi.ProjectionOnX = &PXtest ;
  
  /* Call the callback */
  double x[10], F[10], PX[10];
  int i, n=10;
  for (i =0; i< n ; i++)
  {
    x[i] = i-5;
  }
  vi.F(&vi,n,x,F);
  vi.ProjectionOnX(&vi,x,PX);
  for (i =0; i< n ; i++)
  {
    printf("x[%i]=%f\t",i,x[i]);    printf("F[%i]=%f\t",i,F[i]);    printf("PX[%i]=%f\n",i,PX[i]);
  }
  

    
}
