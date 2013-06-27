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

/*! \file SlidingReducedOrderObserver.hpp
  \brief Discrete-time Sliding Observer
  */

#ifndef SlidingReducedOrderObserver_H
#define SlidingReducedOrderObserver_H

#include "Observer.hpp"
#include "SiconosAlgebraTypeDef.hpp"

class SlidingReducedOrderObserver : public Observer
{
private:
  /** serialization hooks */
  ACCEPT_SERIALIZATION(SlidingReducedOrderObserver);

  /** default constructor */
  SlidingReducedOrderObserver() {};

protected:

  /** the vector defining the measurements (\f$ y = Cx \f$) */
  SP::SiconosMatrix _C;

  /** matrix multiplying the innovation term */
  SP::SiconosMatrix _L;

  double _theta;

  // clumsy hack to do nothing the first time this Observer is called
  bool _pass;

public:

  /** Constructor with the standard interface
   * \param sensor the SP::ControlSensor that feed us with measurements
   * \param xHat0 the initial guess for the state
   */
  SlidingReducedOrderObserver(SP::ControlSensor sensor, const SiconosVector& xHat0):
    Observer(SLIDING_REDUCED_ORDER, sensor, xHat0) {}

  /** Constructor with all the data
   * \param sensor the sensor that feeds the Observer
   * \param xHat0 the initial guess for the state
   * \param C observation matrix
   * \param L gain matrix
   */
  SlidingReducedOrderObserver(SP::ControlSensor sensor, const SiconosVector& xHat0, SP::SiconosMatrix C, SP::SiconosMatrix L):
    Observer(SLIDING_REDUCED_ORDER, sensor, xHat0), _C(C), _L(L) {}

  /** Update the estimate at each event
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
DEFINE_SPTR(SlidingReducedOrderObserver)
#endif
