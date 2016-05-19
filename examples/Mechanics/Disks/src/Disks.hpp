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
