/* Siconos-sample version 3.1.0, Copyright INRIA 2005-2009.
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
 *
 */


#ifndef Disks_h
#define Disks_h

// 2D
#define NDOF 3

// WALLS, TOP and GROUND
#define WALL 100
#define TOP 100
#define GROUND 0

// DEFAULT PLANS : a ground and two walls to support a crystal

// CRYSTAL SIZE
#ifndef Ll
#define Ll 256
#endif

#define Rr 1

#define COSPI6  0.866025403784439
#define SINPI6  0.5
#define TANPI6  0.577350269189626 // tan(pi/6)

#define SY 3.73205080756888  // ((cos(a)+1)/(cos(a)*sin(a)) - tan(a)) a=pi/6, R=1
#define SYL 1/TANPI6


// Plan1
#define P1A COSPI6
#define P1B -SINPI6
#define P1C (SY+(Ll-1)*SYL-Rr)*P1B

// Plan2
#define P2A COSPI6
#define P2B SINPI6
#define P2C (SY+(Ll-1)*SYL-Rr)*P2B


#define GROUND_ID -1
#define MAX_RADIUS INFINITY

#include "SiconosBodies.hpp"

class Disks : public SiconosBodies, public std11::enable_shared_from_this<Disks>
{
  void init();
  void compute();
};

TYPEDEF_SPTR(Disks);

#endif //Disks_h
