#include "g2c.h"
#include <stdio.h>



#define F77NAME(x) x##_

typedef void (*ptr_fex)(integer *, doublereal *, doublereal *, doublereal *);
typedef void (*ptr_jex)(integer *, doublereal *, doublereal *, integer* , integer *,  doublereal *, integer *);


extern "C" void F77NAME(f1)(integer *sizeOfX, doublereal *time, doublereal *x, doublereal *xdot);
extern "C" void F77NAME(jac1)(integer *sizeOfX, doublereal *time, doublereal *x, integer* ml, integer *mu,  doublereal *jacob, integer *nrowpd);

extern "C" void F77NAME(dlsode)(ptr_fex, integer * neq, doublereal * y, doublereal *t, doublereal *tout, integer * itol, doublereal * rtol, doublereal *atol, integer * itask, integer *istate, integer * iopt, doublereal * rwork, integer * lrw, integer * iwork, integer * liw, ptr_jex, integer * mf);


int main(void)
{

  const doublereal tout1 = 1.39283880203;
  const doublereal dtout = 2.214773875;
  integer nerr = 0;
  integer  itol = 1;
  doublereal  rtol = 0.0;
  doublereal  atol = 1.0e-6;
  integer  lrw = 697;
  integer  liw = 45;
  integer  iopt = 0;
  //
  //First problem
  //
  integer neq = 2;
  integer nout = 4;


  doublereal t = 0.0;
  doublereal y[2] = {2.0, 0.0};
  doublereal  rwork[697];
  integer  iwork[45];

  integer   itask = 1;
  integer  istate = 1;
  doublereal  tout = tout1;
  doublereal  ero = 0.0e0;

  integer mf ;
  doublereal hu;
  integer nqu;

  for (int meth = 1; meth <= 2; meth++)
  {
    for (int miter = 0; miter <= 3; miter++)
    {
      mf = 10 * meth + miter;
      printf(" Solution with mf = %3d\n", mf);
      printf("     t               x               xdot       nq      h\n");

      t = 0.0;
      y[0] = 2.0;
      y[1] = 0.0;


      itask = 1;
      istate = 1;
      tout = tout1;
      ero = 0.0e0;

      for (int i = 1; i <= nout; i++)
      {
        F77NAME(dlsode)(F77NAME(f1), &neq, y, &t, &tout, &itol, &rtol, &atol, &itask, &istate, &iopt, rwork, &lrw, iwork, &liw, F77NAME(jac1), &mf);
        hu = rwork[10];
        nqu = iwork[13];
        printf("%15.4e\t%16.4e\t%14.2e\t%5d\t%14.2e\n", t, y[0], y[1], nqu, hu);
        tout = tout + dtout;
      }

    }
  }


}
