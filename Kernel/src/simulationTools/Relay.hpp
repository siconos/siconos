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
/*! \file Relay.hpp
  \brief Linear Complementarity Problem formulation and solving
*/

#ifndef Relay_H
#define Relay_H

#include "LinearOSNS.hpp"
TYPEDEF_SPTR(RelayProblem)

class Simulation;
/** Formalization and Resolution of a Linear Complementarity Problem (Relay)
 
   \author SICONOS Development Team - copyright INRIA
   \version 3.0.0.
   \date (Creation) Apr 26, 2004
 
  \section Relayintro Aim of the Relay class
 
  This class is devoted to the formalization and the resolution of the
  Relay NonSmooth problems.
  \f[
  w =  q + M z
  \f]
  \f[
  w \geq 0, z \geq 0,  z^{T} w =0
  \f]
  where
     - \f$ w \in R^{n} \f$  and \f$z \in R^{n} \f$ are the unknowns,
     - \f$ M \in R^{n \times n } \f$  and \f$q \in R^{n} \f$
      
  \todo : add "recover" function to start from old values of z and w.
  \todo : review this introduction ...
*/
class Relay : public LinearOSNS
{

protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(Relay);


  /** contains the vector lb (lower bounds) of a Relay system */
  SP::SiconosVector _lb;

  /** contains the vector ub (upper bounds) of a Relay system */
  SP::SiconosVector _ub;

  /** contains the numerics proble for Relay system */
  SP::RelayProblem _numerics_problem;

  /** nslaw effects : visitors experimentation
   */

  struct _BoundsNSLEffect;
  friend struct _BoundsNSLEffect;



public:

  /** xml constructor
   *  \param SP::OneStepNSProblemXML : the XML linked-object
   */
  Relay(SP::OneStepNSProblemXML onestepnspbxml):
    LinearOSNS(onestepnspbxml) {};

  /** constructor from data
   *  \param int id of numerics solver
   */
  Relay( int newNumericsSolverId = SICONOS_RELAY_LEMKE);

  /** destructor
   */
  ~Relay();
  // --- lb ---
  /** get the value of lb, the   lower bounds of the Relay system
   *  \return SiconosVector
   *  \warning: SiconosVector is an abstract class => can not be an
   *  lvalue => return SiconosVector
   */
  inline const SiconosVector getLb() const
  {
    return *_lb;
  }

  /** get lb, the lower bounds of the Relay system
   *  \return pointer on a SiconosVector
   */
  inline SP::SiconosVector lb() const
  {
    return _lb;
  }

  /** set lb to pointer newPtr
   *  \param SP::SiconosVector  newPtr
   */
  inline void setLb(SP::SiconosVector newPtr)
  {
    _lb = newPtr;
  }


  // --- ub ---
  /** get the value of ub, the  upper bounds of the Relay system
   *  \return SiconosVector
   *  \warning: SiconosVector is an abstract class => can not be an
   *  lvalue => return SiconosVector
   */
  inline const SiconosVector getUb() const
  {
    return *_ub;
  }

  /** get lb, the lower bounds of the Relay system
   *  \return pointer on a SiconosVector
   */
  inline SP::SiconosVector ub() const
  {
    return _ub;
  }

  /** set lb to pointer newPtr
   *  \param SP::SiconosVector  newPtr
   */
  inline void setUb(SP::SiconosVector newPtr)
  {
    _ub = newPtr;
  }

  void initialize(SP::Simulation sim) ;

  /** Compute the unknown z and w and update the Interaction (y and lambda )
   *  \param double : current time
   *  \return int, information about the solver convergence.
   */
  int compute(double);

  /** print the data to the screen
   */
  void display() const;

  /** encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param OneStepNSProblem* : the one step problem which must be converted
   * \return a pointer on the problem if it is of the right type, NULL otherwise
   */
  static Relay* convert(OneStepNSProblem* osnsp);

};

#endif // Relay_H
