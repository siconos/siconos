#include <math.h>

#include <math.h>
void Inertia(M, q)
double M[4];
double q[2];
{
  double t2;
  {
    t2 = cos(q[1]);
    M[0] = 3.0 + 2.0 * t2;
    M[1] = t2 + 1.0;
    M[2] = M[1];
    M[3] = 1.0;
    return;
  }
}

