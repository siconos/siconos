/* Siconos-Kernel, Copyright INRIA 2005-2015
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

#include <cassert>
#include <iostream>
#include "NormalConeNSL.hpp"

#include "SimpleMatrix.hpp"
#include "SiconosVector.hpp"


// Default (private)
NormalConeNSL::NormalConeNSL(): NonSmoothLaw(), _H(std11::shared_ptr<SimpleMatrix>()), _K(std11::shared_ptr<SiconosVector>())
{}

NormalConeNSL::NormalConeNSL(unsigned size, SP::SimpleMatrix H, SP::SiconosVector K):
  NonSmoothLaw(size), _H(H), _K(K)
{
assert(H->size(1) == size &&
      "NormalConeNSL::NormalConeNSL - the number of columns in H and the declared size are not equal, check your code !");
}

NormalConeNSL::~NormalConeNSL()
{}

bool NormalConeNSL::isVerified(void) const
{
  bool res = false;
  // to do
  return res;
}

void NormalConeNSL::display() const
{
  std::cout << "------------------------------------" <<std::endl;
  std::cout << "____ data of the NormalConeNSL" <<std::endl;
  std::cout << "| nSLawSize : " << _size <<std::endl;
  std::cout << "| H : " << std::endl;
  _H->display();
  std::cout << "| K : " << std::endl;
  _K->display();
  std::cout << "____________________________" <<std::endl;
  std::cout << "------------------------------------" <<std::endl;
}
