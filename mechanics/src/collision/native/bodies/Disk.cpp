/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2020 INRIA.
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

#include "Disk.hpp"

void Disk::MassSetup()
{
  _mass.reset(new SimpleMatrix(_ndof, _ndof));
  //  mass->resize(ndof,ndof);
  _mass->zero();
  (*_mass)(0, 0) = (*_mass)(1, 1) = massValue;
  (*_mass)(2, 2) = massValue * radius * radius / 2.;
}

Disk::Disk(double r, double m,
           SP::SiconosVector qinit,
           SP::SiconosVector vinit)
  : CircularDS(r, m, qinit, vinit)
{
  MassSetup();
}

Disk::~Disk()
{}
