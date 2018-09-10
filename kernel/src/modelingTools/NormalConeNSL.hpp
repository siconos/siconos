/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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
/*! \file NormalConeNSL.hpp
    \brief formalization of the NormalCone nonsmooth law
*/
#ifndef NORMALCONENSLAW_H
#define NORMALCONENSLAW_H

#include "NonSmoothLaw.hpp"

/** NormalCone NonSmoothLaw
 
  This class formalizes a nonsmooth law in the form of a normal cone inclusion i.e.

  \rst
  
  .. math::
  
      y \in \mathcal{N}_{-P}(-\lambda),
      
  \endrst

  where \f$P\f$ is the polytopic set. This is a generalization of the RelayNSL law,
  where the set \f$P\f$ is a scaled box. Note that there exists an inverse of the
  previous relation in the form

  \rst

  .. math::

      -\lambda \in \partial \sigma_{-P} (y),
   
  \endrst

  with \f$\sigma_{-P}\f$ the support function of \f$-P\f$ and \f$\partial \sigma_{-P}\f$
  the subdifferential of this support function.
 
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
