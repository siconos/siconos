/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/

#include <cmath>
#include "DiskPlanR.hpp"
#include "Disk.hpp"
#include <BlockVector.hpp>

void DiskPlanR::init(double r, double A, double B, double C,
                     double xC, double yC, double w)
{
  sqrA2pB2 = hypot(A, B);

  if (width == 0 || width == std::numeric_limits<double>::infinity())
    finite = false;
  else
    finite = true;
  if (fabs(A * xC + B * yC + C) > std::numeric_limits<double>::epsilon())
  {
    if (A == 0)
      // By+C=0
    {
      yCenter = -C / B;
    }
    else if (B == 0)
      // Ax+C=0
    {
      xCenter = -C / A;
    }
    else
      // Ax+By+C=0
    {
      if (xCenter != 0)
        yCenter = - (A * xCenter + C) / B;
      else
        xCenter = - (B * yCenter + C) / A;
    }
  }

  AC = A * C;
  B2 = B * B;
  A2 = A * A;
  AB = A * B;
  BC = B * C;
  x1 = xCenter + halfWidth * fabs(B) / sqrA2pB2;
  y1 =  yCenter - copysign(1, B) * halfWidth * A / sqrA2pB2;
  x2 = xCenter - halfWidth * fabs(B) / sqrA2pB2;
  y2 = yCenter + copysign(1, B) * halfWidth * A / sqrA2pB2;

  halfWidth = 0.5 * width;

}

DiskPlanR::DiskPlanR(double r, double A, double B, double C) :
  LagrangianScleronomousR(),
  r(r), A(A), B(B), C(C)
{
  init(r, A, B, C, 0, 0, std::numeric_limits<double>::infinity());
}

DiskPlanR::DiskPlanR(double r, double A, double B, double C,
                     double xC, double yC, double w) :
  LagrangianScleronomousR(),
  r(r), A(A), B(B), C(C),
  xCenter(xC), yCenter(yC), width(w)
{
  init(r, A, B, C, xCenter, yCenter, width);
}

DiskPlanR::DiskPlanR(double r, double xa, double ya, double xb, double yb)
{

  double d = xa * yb - xb * ya;
  init(r, (yb - ya) / d, (xb - xa) / d, 1, (xa + xb) / 2, (ya + yb) / 2, hypot(xa - xb, ya - yb));
}

double DiskPlanR::distance(double x, double y, double rad) const
{
  if (finite)
  {
    double x0 = - (AC - B2 * x + AB * y) / (A2 + B2);
    double y0 = - (BC - A2 * y + AB * x) / (A2 + B2);

    if (hypot(xCenter - x0, yCenter - y0) >= halfWidth)
    {
      // ... and no jach...
      double r = fmin(hypot(x - x1, y - y1), hypot(x - x2, y - y2));
      return r - rad;
    }
  }

  return (fabs(A * x + B * y + C) / sqrA2pB2 - rad);

}

/* called compute H, but only the gap function is needed! */
void DiskPlanR::computeh(SiconosVector& q, SiconosVector& z, SiconosVector& y)
{
  double q_0 = q(0);
  double q_1 = q(1);

  y(0) = distance(q_0, q_1, r);

}

/* Note :

   Only the first component (normal part) of H is meaningfull for siconos
   see comment in Interaction.cpp/getYRef

   The normal part (gap function) is used in updateIndexSet.

   The second part is the same as the one of the linear mapping (i.e
   what would be done with a LagrangianLinearTIR)

    t
   q  = [ x y theta ]
        t
   y = H  q + b
            2     2
   c = sqrt(A   + B  )

    t
   H =|  A/c    B/c     0|
      | -B/c    A/c    -r|

   b = C / c

   y[0]  = A/c x + B/c y + C/c >= 0

   Only one side.

   Both side with g(q) = | A/c x + B/c y + C/c | as the gap function
                                      t
   Hnl(q) = |   g(q)          |    = H  q   +  |   g(q) - (A/c x + B/c y) |
            | -B/c    A/c  -R |                |   0                      |


  so:


           [           |C + A*x + B*y|            ]
           [           ---------------            ]
           [                _________             ]
           [               /  2    2              ]
           [             \/  A  + B               ]
   H(q) =  [                                      ]
           [               A*y            B*x     ]
           [-r*theta + ------------ - ------------]
           [              _________      _________]
           [             /  2    2      /  2    2 ]
           [           \/  A  + B     \/  A  + B  ]


           [A*sign(C + A*x + B*y)  B*sign(C + A*x + B*y)    ]
           [---------------------  ---------------------  0 ]
           [        _________              _________        ]
           [       /  2    2              /  2    2         ]
           [     \/  A  + B             \/  A  + B          ]
Jach(q) =  [                                                ]
           [         -B                     A               ]
           [    ------------           ------------       -r]
           [       _________              _________         ]
           [      /  2    2              /  2    2          ]
           [    \/  A  + B             \/  A  + B           ]


*/


void DiskPlanR::computeJachq(SiconosVector& q, SiconosVector& z)
{

  SimpleMatrix *g = (SimpleMatrix *) _jachq.get();

  double x = q(0);
  double y = q(1);

  double D1 = A * x + B * y + C;
  double signD1 = copysign(1, D1);

  (*g)(0, 0) = A * signD1 / sqrA2pB2;
  (*g)(1, 0) = -B * signD1 / sqrA2pB2;
  (*g)(0, 1) = B * signD1 / sqrA2pB2;
  (*g)(1, 1) = A * signD1 / sqrA2pB2;
  (*g)(0, 2) = 0;
  (*g)(1, 2) = -r;
}

bool DiskPlanR::equal(double pA, double pB, double pC, double pr) const
{
  return (A == pA && B == pB && C == pC && r == pr);
}

bool DiskPlanR::equal(double pA, double pB, double pC, double pr,
                      double pXc, double pYc, double pw) const
{
  if (finite)
    return (A == pA && B == pB && C == pC && r == pr &&
            pXc == xCenter && pYc == yCenter && pw == width);
  else
    return equal(pA, pB, pC, pr);
}

bool DiskPlanR::equal(const DiskPlanR& odpr) const
{
  if (finite)
    return (equal(odpr.getA(), odpr.getB(), odpr.getC(),
                  odpr.getRadius()));
  else
    return (equal(odpr.getA(), odpr.getB(), odpr.getC(),
                  odpr.getRadius(), odpr.getXCenter(),
                  odpr.getYCenter(), odpr.getWidth()));
}
