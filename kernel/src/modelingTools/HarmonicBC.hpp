/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
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
#ifndef HARMONICBC_HPP
#define HARMONICBC_HPP

#include "BoundaryCondition.hpp"
#include "SiconosVector.hpp"

namespace siconos::algebra {
}

namespace siconos::modeling {
/**\brief This class models a simple harmonic boundary conditions for
 *   prescribing the velocities in a Dynamical System. A simple
 *   boundary condition is considered to fix a component \f$ j \f$ of
 *   the velocity vector, i.e., \f$ v_j(t) = a +  b cos( \omega t+ \phi) \f$.
 *
 */
class HarmonicBC : public BoundaryCondition {
 public:
  /** Constructor
   * \param newVelocityIndices the indices of the velocity subjected to prescribed velocities
   * \param a constant value for additive term of the prescribed velocity
   * \param b constant value for multiplicative term of the prescribed velocity
   * \param omega frequency
   * \param phi phase
   */
  HarmonicBC(Indices&& newVelocityIndices, double a, double b, double omega, double phi)
    : BoundaryCondition(std::move(newVelocityIndices)), _a(a), _b(b), _omega(omega), _phi(phi){};

  HarmonicBC(Indices&& newVelocityIndices, std::shared_ptr<siconos::algebra::SiconosVector> a,
             std::shared_ptr<siconos::algebra::SiconosVector> b,
             std::shared_ptr<siconos::algebra::SiconosVector> omega,
             std::shared_ptr<siconos::algebra::SiconosVector> phi);

  /** destructor */
  virtual ~HarmonicBC() noexcept = default;

  /** default function to compute the precribed velocities
   *
   *  \param  time : the current time
   */
  virtual void computePrescribedVelocity(double time);

 protected:
  ACCEPT_SERIALIZATION(HarmonicBC);

  /** Constant additive term of the prescribed velocity  */
  double _a = 0.;
  /** Constant multiplicative term of the prescribed velocity  */
  double _b = 0;
  ;
  /** Constant frequency  */
  double _omega = 0;
  ;
  /** Constant phase  */
  double _phi = 0;
  ;

  /** Constant additive term of the prescribed velocity  */
  std::shared_ptr<siconos::algebra::SiconosVector> _aV{nullptr};
  /** Constant multiplicative term of the prescribed velocity  */
  std::shared_ptr<siconos::algebra::SiconosVector> _bV{nullptr};
  /** Constant frequency  */
  std::shared_ptr<siconos::algebra::SiconosVector> _omegaV{nullptr};
  /** Constant phase  */
  std::shared_ptr<siconos::algebra::SiconosVector> _phiV{nullptr};

  //   /*Link to the precribed DynamicalSystem*/
  //   std::shared_ptr<DynamicalSystem> _DS;
};
}  // namespace siconos::modeling
#endif  // HARMONICBC_HPP
