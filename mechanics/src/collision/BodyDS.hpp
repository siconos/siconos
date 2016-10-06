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

/*! \file BodyDS.hpp
  \brief Definition of an abstract body
*/


#ifndef BodyDS_h
#define BodyDS_h

#include "MechanicsFwd.hpp"
#include "NewtonEulerDS.hpp"
#include <Question.hpp>

class BodyDS : public NewtonEulerDS, public std11::enable_shared_from_this<BodyDS>
{
protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(BodyDS);

  SP::SiconosContactor _contactor;

public:

  BodyDS(SP::SiconosVector position,
         SP::SiconosVector velocity,
         double mass);

  virtual ~BodyDS();

  void setContactor(SP::SiconosContactor contactor)
    { _contactor = contactor; }

  SP::SiconosContactor contactor() const { return _contactor; }

  /** visitors hook
   */
  ACCEPT_BASE_STD_VISITORS(NewtonEulerDS);
};

#endif /* BodyDS_h */
