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

/*! \file SiconosBodies.hpp
  \brief SiconosBodies class - model + plans + space filter
*/
#ifndef SiconosBodies_hpp
#define SiconosBodies_hpp

#include "Model.hpp"
#include "SpaceFilter.hpp"

/** SiconosBodies : a Siconos Model, some plans and space filtering capabilities
 */

class SiconosBodies
{

protected:

  SP::SiconosMatrix plans_;
  SP::Model model_;
  SP::SpaceFilter playground_;

public:

  virtual void init() = 0;

  virtual void compute();

  SP::Model model()
  {
    return model_;
  }

  SP::SiconosMatrix getPlansPtr()
  {
    return plans_;
  }

};

TYPEDEF_SPTR(SiconosBodies);

#endif // SiconosBodies_hpp
