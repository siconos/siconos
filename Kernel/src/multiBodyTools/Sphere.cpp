/* Siconos-Example version 3.0.0, Copyright INRIA 2005-2008.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth FLOOR, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
 *
 */

#include "Sphere.h"

using namespace std;

Sphere::Sphere(double R, double m, const SiconosVector& q0, const SiconosVector& v0): LagrangianDS(q0, v0), Radius(R), mass(m), nDof(6)
{
  // Gravity
  double g = 9.81;
  // Mass matrix (the same for all beads)
  SimpleMatrix Mass(nDof, nDof);
  Mass(0, 0) = Mass(1, 1) = Mass(2, 2) = m;    ;
  Mass(3, 3) = Mass(4, 4) = Mass(5, 5) = 3. / 5 * m * Radius * Radius;
  // External forces = gravity
  SimpleVector weight(nDof);
  weight(2) = -mass * g;

  setFExt(weight);
  setMass(Mass);
}

Sphere::~Sphere()
{}

void Sphere::draw()
{
  //   float theta;//,phi,psi;
  //   float pos[3];
  //   float R[12];
  //   // Translation
  //   pos[0]=(*q[0])(0);
  //   pos[1]=(*q[0])(1);
  //   pos[2]=(*q[0])(2);

  //   // Rotation
  //   theta=(*q[0])(3);
  // //   phi=getQ()(4);
  // //   psi=getQ()(5);

  //   R[0] = 1.0; R[1] = 0.0; R[2] = 0.0;
  //   R[4] = 0.0; R[5] = cos(theta); R[6] = -sin(theta);
  //   R[8] = 0.0; R[9] = sin(theta); R[10] = cos(theta);
  //   R[3]=R[7]=R[11]=0.0;

  //   dsSetTexture (DS_WOOD);
  //   dsSetColor (1,0.8f,0.6f);

  //   dsDrawSphere (pos,R, Radius);

  // static void drawSphere(GLdouble radius, double x, double y, double z)
  // {
}
