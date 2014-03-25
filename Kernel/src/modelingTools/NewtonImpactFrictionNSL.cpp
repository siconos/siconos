/* Siconos-Kernel, Copyright INRIA 2005-2012.
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
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
*/
#include "NewtonImpactFrictionNSL.hpp"

#include <iostream>

// Default (private)
NewtonImpactFrictionNSL::NewtonImpactFrictionNSL():
  NonSmoothLaw(), _en(0.0), _et(0.0), _mu(0.0)
{}

NewtonImpactFrictionNSL::NewtonImpactFrictionNSL(unsigned int size):
  NonSmoothLaw(size), _en(0.0), _et(0.0), _mu(0.0)
{}

NewtonImpactFrictionNSL::NewtonImpactFrictionNSL(double newEn, double newEt, double newMu, unsigned int newSize):
  NonSmoothLaw(newSize), _en(newEn), _et(newEt), _mu(newMu)
{}

NewtonImpactFrictionNSL::~NewtonImpactFrictionNSL()
{}

bool NewtonImpactFrictionNSL::isVerified() const
{
  bool res = false;
  // to do
  RuntimeException::selfThrow("NewtonImpactFrictionNSL:: isVerified, not yet implemented!");
  return res;
}

void NewtonImpactFrictionNSL::display() const
{
  std::cout << "=== Newton impact-friction non-smooth law data display ===" <<std::endl;
  std::cout << " Normal Newton coefficient of restitution: " << _en <<std::endl;
  std::cout << " Tangential Newton coefficient of restitution: " << _et <<std::endl;
  std::cout << "Friction coefficient: " << _mu <<std::endl;
  std::cout << "==========================================================" <<std::endl;
}
