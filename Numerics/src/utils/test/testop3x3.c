
/* random is POSIX but not C99...*/
#define _BSD_SOURCE

#include <stdlib.h>
#include <stdio.h>
#include <LA.h>
#include <assert.h>

#define WITH_TIMER

#include "timers.h"

#include "op3x3.h"



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


  DECL_TIMER_TICK(t1);
  DECL_TIMER_TICK(t2);
  DECL_TIMER_TICK(t3);
  DECL_TIMER_TICK(t4);
  DECL_TIMER_TICK(t5);
  DECL_TIMER_TICK(t6);
  DECL_TIMER_TICK(t7);
  DECL_TIMER_TICK(t8);

  unsigned int i, k, k1, k2;

  DECL_TIMER(t);


  /* b += a */
  cpy3x3(ia, a);
  cpy3x3(ib, b);

  START_TIMER(t);
  for (i = 0; i < 1000000; ++i)
    DAXPY(9, 1., a, 1, b, 1);
  STOP_TIMER(t);
  GET_ELAPSED(t, t1);
  PRINT_ELAPSED(t);

  cpy3x3(b, c);

  cpy3x3(ia, a);
  cpy3x3(ib, b);
  START_TIMER(t);
  for (i = 0; i < 1000000; ++i)
    add3x3(a, b);
  STOP_TIMER(t);
  GET_ELAPSED(t, t2);
  PRINT_ELAPSED(t);

  assert(equal3x3(b, c));


#ifdef WITH_TIMER
  printf("add3x3/DAXPY:%g\n", t2 / t1);
  assert(t2 < t1);
#endif

  /* c += a*b */
  cpy3x3(ia, a);
  cpy3x3(ib, b);
  cpy3x3(ib, c);
  START_TIMER(t);
  for (i = 0; i < 1000000; ++i)
  {
    DGEMM(LA_NOTRANS, LA_NOTRANS, 3, 3, 3, 1, a, 3, b, 3, 1, c, 3);
  }
  STOP_TIMER(t);
  GET_ELAPSED(t, t3);

  print3x3(c);
  cpy3x3(c, d);


  cpy3x3(ia, a);
  cpy3x3(ib, b);
  cpy3x3(ib, c);
  START_TIMER(t);
  for (i = 0; i < 1000000; ++i)
  {
    mmp3x3(a, b, c);
  }
  STOP_TIMER(t);
  GET_ELAPSED(t, t4);

  assert(equal3x3(c, d));

  PRINT_ELAPSED(t);

#ifdef WITH_TIMER
  printf("mmp3x3/DGEMM:%g\n", t4 / t3);
  assert(t4 < t3);
#endif

  print3x3(c);



  cpy3x3(ia, a);
  cpy3x3(ib, b);
  cpy3x3(ib, c);

  START_TIMER(t);
  for (i = 0; i < 1000000; ++i)
  {
    r = det3x3(a);
  }
  STOP_TIMER(t);

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

  int ipiv[3];
  int info;

  START_TIMER(t);
  for (i = 0; i < 1; ++i)
  {
    cpy3x3(v, x);
    DGESV(3, 1, a, 3, ipiv, x, 3, info);
    v[0] += 0.0001;
  }
  STOP_TIMER(t);
  GET_ELAPSED(t, t5);
  PRINT_ELAPSED(t);
  printf("sol DGESV:\n");
  print3(x);
  cpy3x3(x, d);




  cpy3x3(ia, a);
  cpy3(iv, v);

  START_TIMER(t);
  for (i = 0; i < 1; ++i)
  {
    cpy3x3(v, x);
    solv3x3(a, x, v);
    v[0] += 0.0001;
  }
  STOP_TIMER(t);
  GET_ELAPSED(t, t6);

  printf("sol solv3x3:\n");
  print3(x);

  sub3x3(x, d);

  assert(hypot3(d) <= 1e-5);

  mv3x3(a, x, c);

  printf("a*x:\n");
  print3(c);


  sub3x3(iv, c);
  printf("v - a*x:\n");
  print3(c);

  assert(hypot3(c) <= 1e-5);


  PRINT_ELAPSED(t);

#ifdef WITH_TIMER
  printf("solv3x3/DGESV:%g\n", t6 / t5);
  assert(t6 < t5);
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

  START_TIMER(t);
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

    DGESV(3, 1, a, 3, ipiv, x, 3, info);

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

      sub3x3(B + 3 * i, c);

      /* printf("hypot3(c)=%g\n",hypot3(c));*/

      assert(hypot3(c) <= 1e-5);
    }
  }
  STOP_TIMER(t);
  GET_ELAPSED(t, t7);
  k1 = k;

  k = 0;
  START_TIMER(t);
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
    solv3x3(a, x, B + 3 * i);

    /*
    printf("x:\n");
    print3(x);
    */

    if (! isnan3(x))
    {
      k++;

      mv3x3(A + 9 * i, x, c);

      /*
      printf("Ax:\n");

      print3(c);
      */

      sub3x3(B + 3 * i, c);


      /*printf("hypot3(c)=%g\n",hypot3(c));*/

      assert(hypot3(c) <= 1e-5);
    }
  }
  STOP_TIMER(t);
  GET_ELAPSED(t, t8);
  k2 = k;


#ifdef WITH_TIMER
  printf("DGESV/solv3x3 (rand):%g\n", t8 / t7);
  printf("Number of pb DGESV and solv3x3 (rand):%d,%d\n", k1, k2);
  assert(t8 < t7);
#endif


  return 0;

}
