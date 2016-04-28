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
/*! \file NormalConeNSL.hpp
    \brief formalization of the NormalCone nonsmooth law
*/
#ifndef NORMALCONENSLAW_H
#define NORMALCONENSLAW_H

#include "NonSmoothLaw.hpp"

/** NormalCone NonSmoothLaw
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.8.0.
 *  \date Feb 10, 2015
 *
 * This class formalizes a nonsmooth law in the form of a normal cone inclusion i.e.
 * \f[
 * y \in \mathcal{N}_{-P}(-\lambda),
 * \f]
 * where \f$P\f$ is the polytopic set. This is a generalization of the RelayNSL law,
 * where the set \f$P\f$ is a scaled box. Note that there exists an inverse of the
 * previous relation in the form
 * \f[
 * -\lambda \in \partial \sigma_{-P} (y),
 *  \f]
 *  with \f$\sigma_{-P}\f$ the support function of \f$-P\f$ and \f$\partial \sigma_{-P}\f$
 *  the subdifferential of this support function.
 *
 */
class NormalConeNSL : public NonSmoothLaw
{

private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(NormalConeNSL);

  /** matrix in the (H-K)-representation of the polytope */
  SP::SimpleMatrix _H;

  /** vector in the (H-K)-representation of polytope */
  SP::SiconosVector _K;

  /** default constructor
   */
  NormalConeNSL();

public:

  /** constructor with the value of the NormalConeNSL attributes
  *  \param size size of the NonSmoothLaw
  *  \param H matrix in the (H-K)-representation of the polytope P
  *  \param K vector in the (H-K)-representation of the polytope P
  */
  NormalConeNSL(unsigned size, SP::SimpleMatrix H, SP::SiconosVector K);

  virtual ~NormalConeNSL();

  /** get H
   * \return a reference to the H matrix
   */
  inline SimpleMatrix& H() { return *_H; };

  /** get K
   * \return a reference to the K vector
   */
  inline SiconosVector& K() { return *_K; };

  /** check the ns law to see if it is verified
  *  \return true if the NS Law is verified, false otherwise
  */
  virtual bool isVerified() const;

  /** print the data to the screen */
  virtual void display() const;

  /** visitors hook
   */
  ACCEPT_STD_VISITORS();

};

#endif // NORMALCONENSLAW_H
