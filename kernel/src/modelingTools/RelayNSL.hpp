/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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
/*! \file RelayNSL.hpp
    \brief formalization of the Relay nonsmooth law
*/
#ifndef RELAYNSLAW_H
#define RELAYNSLAW_H

#include "NonSmoothLaw.hpp"

/** 
    Relay NonSmoothLaw

    This class formalizes the Relay nonsmooth law  i.e.

    \f[
    -y \in \mathcal{N}_{[lb,ub]}(\lambda),
    \f]
    
    where \f$ lb \f$ is the lower bound and   \f$ ub \f$ is the upper bound of the
    Relay law.
    
    In this default case, the lower bound is set to \f$ lb=-1 \f$ and the upper bound
    ub is set to \f$ ub=1 \f$. We get  the well-known form of the RelayNSL as the
    multivalued sign function, i.e.
    
    \f[
    y \in -\mathcal{N}_{[-1,1]}(\lambda) \Longleftrightarrow \lambda \in -\mbox{sgn} (y)  
    \f]

    where the multi-valued sign function is defined as
    
    \f[

    \mbox{sgn} (y) =
    \left\{ \begin{array}{lcl}
    1 && y >0 \\
    [-1,1] && y =0 \\
    -1 && y <0
    \end{array}\right.

    \f]

    \todo Build the Sgn NonSmoothLaw as the default instance of Relay

*/
class RelayNSL : public NonSmoothLaw {

private:
  ACCEPT_SERIALIZATION(RelayNSL);

  /** represent the lower bound of the Relay */
  double _lb;

  /** represent the upper bound of the Relay*/
  double _ub;

  /** default constructor
   */
  RelayNSL();

public:
  /** constructor with the value of the RelayNSL attributes
   *
   *  \param size size of the NonSmoothLaw
   *  \param lb lower endpoint of the interval, default value is -1.0
   *  \param ub upper endpoint of the interval, default value is 1.0
   */
  RelayNSL(unsigned int size, double lb = -1.0, double ub = 1.0);

  ~RelayNSL();

  /** check the ns law to see if it is verified
   *
   *  \return true if the NS Law is verified, false otherwise
   */
  bool isVerified() const override;

  /** to get lb
   *
   *  \return the value of lb
   */
  inline double lb() const { return _lb; };

  /** to set the lower bound
   *
   *  \param lb the new lower bound
   */
  inline void setLb(double lb) { _lb = lb; };

  /** to get ub
   *
   *  \return the value of ub
   */
  inline double ub() const { return _ub; };

  /** to set ub
   *
   *  \param ub the new upper bound
   */
  inline void setUb(double ub) { _ub = ub; };

  /** print the data to the screen
   */
  void display() const override;

  ACCEPT_STD_VISITORS();
};

#endif // RELAYNSLAW_H
