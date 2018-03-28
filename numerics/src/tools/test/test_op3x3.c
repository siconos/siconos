
/* random is POSIX but not C99...
 * we don't define _BSD_SOURCE but rather _XOPEN_SOURCE */
#define _XOPEN_SOURCE 700

#include <stdlib.h>
#include <stdio.h>
#include "SiconosLapack.h"
#include <assert.h>

#define TIMER_FFTW_CYCLE

#include "timers_interf.h"

#include "op3x3.h"


#ifdef _WIN32
#define random rand
#define srandom srand
#endif


int main()
{
  double ia[9] = { 1., 2., 3.,
                   0., 1., 0.,
                   1., 2., 1.
                 };


  double ib[9] = { 1., 2., 0.,
                   0., 1., 0.,
                   0., 2., 1.
                 };

  double a[9], b[9], c[9], d[9], x[9];
  double iv[3] = { .1, .1, .1 };
  double v[3];
  double r;


  DECL_TIMER(t1);
  DECL_TIMER(t2);
  DECL_TIMER(t3);
  DECL_TIMER(t4);
  DECL_TIMER(t5);
  DECL_TIMER(t6);
  DECL_TIMER(t7);
  DECL_TIMER(t8);
  DECL_TIMER(t9);

  unsigned int i, k;

  /* b += a */
  cpy3x3(ia, a);
  assert(equal3x3(ia, a));


  cpy3x3(ib, b);
  assert(equal3x3(ib, b));

  sub3x3(ia, ib);
  add3x3(ia, ib);
  assert(equal3x3(ia, a));

  START_TIMER(t1);
  for (i = 0; i < 1000000; ++i)
    cblas_daxpy(9, 1., a, 1, b, 1);
//    DAXPY(9, 1., a, 1, b, 1);
  STOP_TIMER(t1);
  PRINT_ELAPSED(t1);

  cpy3x3(b, c);

  cpy3x3(ia, a);
  cpy3x3(ib, b);
  START_TIMER(t2);
  for (i = 0; i < 1000000; ++i)
    add3x3(a, b);
  STOP_TIMER(t2);
  PRINT_ELAPSED(t2);

  assert(equal3x3(b, c));


#ifdef WITH_TIMERS
  printf("add3x3/DAXPY:%g\n", ELAPSED(t2) / ELAPSED(t1));
  assert(ELAPSED(t2) < ELAPSED(t1));
#endif

  /* c += a*b */
  cpy3x3(ia, a);
  cpy3x3(ib, b);
  cpy3x3(ib, c);
  START_TIMER(t3);
  for (i = 0; i < 1000000; ++i)
  {
    //DGEMM(CblasNoTrans, CblasNoTrans, 3, 3, 3, 1, a, 3, b, 3, 1, c, 3);
    cblas_dgemm (CblasColMajor, CblasNoTrans, CblasNoTrans, 3, 3, 3, 1.0, a, 3, b, 3, 1.0, c, 3);
     
  }
  STOP_TIMER(t3);
  PRINT_ELAPSED(t3);

  print3x3(c);
  cpy3x3(c, d);


  cpy3x3(ia, a);
  cpy3x3(ib, b);
  cpy3x3(ib, c);
  START_TIMER(t4);
  for (i = 0; i < 1000000; ++i)
  {
    mmp3x3(a, b, c);
  }
  STOP_TIMER(t4);

  assert(equal3x3(c, d));

  PRINT_ELAPSED(t4);

#ifdef WITH_TIMERS
  printf("mmp3x3/DGEMM:%g\n", ELAPSED(t4) / ELAPSED(t3));
  assert(ELAPSED(t4) < ELAPSED(t3));
#endif

  print3x3(c);



  cpy3x3(ia, a);
  cpy3x3(ib, b);
  cpy3x3(ib, c);

  START_TIMER(t5);
  for (i = 0; i < 1000000; ++i)
  {
    r = det3x3(a);
  }
  STOP_TIMER(t5);

  printf("r:%g\n", r);




  /*
  a=matrix([[1, 0, 1],
           [2, 1, 2],
           [3, 0, 1]])

  v=matrix([[ 0.1],
            [ 0.1],
            [ 0.1]])

  linalg.solve(a,v)
  matrix([[  9.25185854e-18],
          [ -1.00000000e-01],
          [  1.00000000e-01]])
  */
  cpy3x3(ia, a);
  cpy3(iv, v);
  assert(equal3(iv, v));
  sub3(iv, v);
  add3(iv, v);
  assert(equal3(iv, v));

  int ipiv[3];
  int info = 0;

  START_TIMER(t6);
  for (i = 0; i < 1; ++i)
  {
    cpy3(v, x);
    DGESV(3, 1, a, 3, ipiv, x, 3, &info);
    printf("info is ...%i \n ", info);
    v[0] += 0.0001;
  }
  STOP_TIMER(t6);
  PRINT_ELAPSED(t6);
  printf("sol DGESV:\n");
  print3(x);
  cpy3x3(x, d);




  cpy3x3(ia, a);
  cpy3(iv, v);
  START_TIMER(t7);
  for (i = 0; i < 1; ++i)
  {
    cpy3(v, x);
    info = solv3x3(a, x, v);
    v[0] += 0.0001;
  }
  STOP_TIMER(t7);

  printf("v:\n");
  print3(v);

  printf("sol solv3x3:\n");
  print3(x);

  sub3x3(x, d);

  assert(hypot3(d) <= 1e-5);

  mv3x3(a, x, c);

  printf("a*x:\n");
  print3(c);


  sub3(iv, c);
  printf("v - a*x:\n");
  print3(c);

  assert(hypot3(c) <= 1e-5);


#ifdef WITH_TIMERS
  printf("solv3x3/DGESV:%g\n", ELAPSED(t7) / ELAPSED(t6));
  assert(ELAPSED(t7) < ELAPSED(t6));
#endif


  srandom(1);
  double A[9000];
  double* pA = A;
  for (i = 0; i < 1000; ++i)
  {
    for (k = 0; k < 9; k++)
    {
      *pA++ = (double) random() / RAND_MAX;
    }
  }

  double B[3000];
  double* pB = B;
  for (i = 0; i < 1000; ++i)
  {
    for (k = 0; k < 3; k++)
    {
      *pB++ = (double) random() / RAND_MAX;
    }
  }

  START_TIMER(t8);
  k = 0;

  for (i = 0; i < 1000; ++i)
  {

    /*
    printf("A:\n");

    print3x3(A+9*i);
    printf("B:\n");

    print3(B+3*i);
    */

    cpy3(B + 3 * i, x);

    cpy3x3(A + 9 * i, a);
    
    DGESV(3, 1, a, 3, ipiv, x, 3, &info);

    if (! info)
    {
      k++;

      /*
      printf("A:\n");
      print3x3(A+9*i);


      printf("x:\n");
      print3(x);
      */

      mv3x3(A + 9 * i, x, c);
      /*
      printf("Ax:\n");
      print3(c);
      */

      sub3(B + 3 * i, c);

      /* printf("hypot3(c)=%g\n",hypot3(c));*/

      assert(hypot3(c) <= 1e-5);
    }
  }
  STOP_TIMER(t8);
#ifdef WITH_TIMERS
  unsigned int k1;
  k1 = k;
#endif
  k = 0;
  START_TIMER(t9);
  for (i = 0; i < 1000; ++i)
  {

    /*
    printf("A:\n");

    print3x3(A+9*i);
    printf("B:\n");

    print3(B+3*i);
    */

    cpy3x3(A + 9 * i, a);

    cpy3(B + 3 * i, x);
    info = solv3x3(a, x, B + 3 * i);

    /*
    printf("x:\n");
    print3(x);
    */

    //    if (! isnan3(x))
    {
      k++;

      mv3x3(A + 9 * i, x, c);

      /*
      printf("Ax:\n");

      print3(c);
      */

      sub3(B + 3 * i, c);


      /*printf("hypot3(c)=%g\n",hypot3(c));*/

      assert(hypot3(c) <= 1e-5);
    }
  }
  STOP_TIMER(t9);

#ifdef WITH_TIMERS
  unsigned int k2;
  k2 = k;
  printf("DGESV/solv3x3 (rand):%g\n", ELAPSED(t9) / ELAPSED(t8));
  printf("Number of pb DGESV and solv3x3 (rand):%d,%d\n", k1, k2);
  assert(ELAPSED(t9) < ELAPSED(t8));
#endif


  return 0;

}
