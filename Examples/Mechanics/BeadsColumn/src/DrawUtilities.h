/* Siconos-Kernel version 3.0.0, Copyright INRIA 2005-2008.
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
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
*/

/*! \file DrawUtilities.h
  \brief Functions used to draw spheres or planes using openGL
  \author F. Perignon
  \date August 2008
*/

#ifdef WithQGLViewer
#ifndef DrawUtilities_H
#define DrawUtilities_H

/** Toolbox for drawing  */
class DrawUtilities
{
public:

  /** draw a sphere of radius R, centered at point of coordinate (x,y,z)
      c is used to choose the color
  */
  static void drawSphere(double radius, double x, double y, double z, double c);

  /** draw the plane z = ground
   */
  static void drawHorizontalPlane(double ground);

};
#endif

#endif
