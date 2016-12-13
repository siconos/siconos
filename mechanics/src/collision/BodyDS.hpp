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

#include <MechanicsFwd.hpp>
#include <NewtonEulerDS.hpp>
#include <SiconosVisitor.hpp>
#include <SiconosContactor.hpp>

class BodyDS : public NewtonEulerDS,
               public std11::enable_shared_from_this<BodyDS>
{
private:
  BodyDS() : NewtonEulerDS() {};

protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(BodyDS);

  SP::SiconosContactorSet _contactors;
  bool _useContactorInertia;

public:

  BodyDS(SP::SiconosVector position,
         SP::SiconosVector velocity,
         double mass,
         SP::SimpleMatrix inertia = SP::SimpleMatrix());

  virtual ~BodyDS();

  void setUseContactorInertia(bool use) { _useContactorInertia = use; }

  bool useContactorInertia() { return _useContactorInertia; }

  /** Access the contactor set associated with this body.
   * \return A SP::SiconosContactorSet */
  SP::SiconosContactorSet contactors() const { return _contactors; }

  /** Provide a set of contactors to the body.
   * \param c A SP::SiconosContactorSet */
  void setContactors(SP::SiconosContactorSet c) { _contactors = c; }

  /** Make the base position of the contactors equal to the DS q vector.
   * \return a SP::SiconosVector */
  virtual SP::SiconosVector base_position() { return q(); }

  /** visitors hook
   */
  ACCEPT_BASE_VISITORS(NewtonEulerDS);
};

#endif /* BodyDS_h */
