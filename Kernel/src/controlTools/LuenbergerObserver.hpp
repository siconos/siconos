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

/*! \file LuenbergerObserver.hpp
  \brief Classical discrete-time Luenberger Observer
  */

#ifndef LuenbergerObserver_H
#define LuenbergerObserver_H

#include "Observer.hpp"
#include "SiconosAlgebraTypeDef.hpp"

#ifndef FirstOrderLinearDS_H
DEFINE_SPTR(FirstOrderLinearDS)
#endif
#ifndef TimeStepping_H
DEFINE_SPTR(TimeStepping)
#endif
#ifndef Relay_H
DEFINE_SPTR(Relay)
#endif
#ifndef EventsManager_H
DEFINE_SPTR(EventsManager)
#endif

class LuenbergerObserver : public Observer
{
private:
  /** serialization hooks */
  ACCEPT_SERIALIZATION(LuenbergerObserver);

  /** default constructor */
  LuenbergerObserver() {};

protected:

  /** the vector defining the measurements (\f$ y = Cx \f$) */
  SP::SiconosMatrix _C;

  /** matrix describing the relation between the control value and sgn(s) */
  SP::SiconosMatrix _L;

  double _theta;

  // clumsy hack to do nothing the first time this Observer is called
  bool _pass;

public:

  /** Constructor with a TimeDiscretisation, a ControlSensor and an initial estimate of the state.
   * \param sensor the SP::ControlSensor that feed us with measurements
   * \param xHat0 the initial guess for the state
   */
  LuenbergerObserver(SP::ControlSensor sensor, const SiconosVector& xHat0):
    Observer(LUENBERGER, sensor, xHat0) {}

  /** Constructor with all the data
   * \param sensor the ControlSensor that feeds the Observer
   * \param xHat0 the initial guess for the state
   * \param C the observation matrix
   * \param L the gain matrix
   */
  LuenbergerObserver(SP::ControlSensor sensor, const SiconosVector& xHat0, SP::SiconosMatrix C, SP::SiconosMatrix L):
    Observer(LUENBERGER, sensor, xHat0), _C(C), _L(L) {}

  /** Compute the new control law at each event
  */
  virtual void process();

  /** Initialization
   * \param m the Model
   */
  virtual void initialize(const Model& m);

  /** Set the L matrix
   * \param L the new L matrix
  */
  void setL(const SiconosMatrix& L);

  /** Set the L matrix
   * \param L the new L matrix
   */
  inline void setLPtr(SP::SiconosMatrix L)
  {
    _L = L;
  };

  /** Set the C matrix
   * \param c the new C matrix
  */
  void setC(const SiconosMatrix& C);

  /** Set the C matrix
   * \param C the new C matrix
   */
  inline void setcPtr(SP::SiconosMatrix C)
  {
    _C = C;
  };

};
DEFINE_SPTR(LuenbergerObserver)
#endif
