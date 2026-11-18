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
/*! \file NewtonEulerVelocityAngularVelocityR.hpp
 *  \brief Contact relation between Newton-Euler bodies: relative velocity and tangential
 *  relative angular velocity in the contact frame (Coulomb friction and rolling resistance).
 */
#ifndef NEWTONEULERVELOCITYANGULARVELOCITYR_HPP
#define NEWTONEULERVELOCITYANGULARVELOCITYR_HPP

#include "NewtonEuler1DR.hpp"

namespace siconos::modeling {

/** Unilateral contact with friction and rolling resistance between body 1 and body 2 (or
  body 1 and a fixed obstacle): relative velocity at the contact point and tangential
  components of the relative angular velocity, in the contact frame.

  The geometry (contact points \f$ P_1, P_2 \f$, normal \f$ n \f$, gap \f$ g \f$), the
  contact frame \f$ R_c = [\, n \;\; t_1 \;\; t_2 \,]^\top \f$ and the matrix \f$ H_c(q) \f$
  are defined in NewtonEuler1DR. With \f$ \omega_i = R_i \Omega_i \f$ the angular
  velocities in the inertial frame, here \f$ m = 5 \f$ and
  \f[
    \dot y = \begin{bmatrix} R_c\, u \\ t_1^\top (\omega_1 - \omega_2) \\
                             t_2^\top (\omega_1 - \omega_2) \end{bmatrix} = H(q)\, v,
    \qquad
    H(q) = \begin{bmatrix}
      R_c & R_c [x_{G_1} - P]_\times R_1 & -R_c & -R_c [x_{G_2} - P]_\times R_2 \\
      0 & (R_c R_1)_{2:3,:} & 0 & -(R_c R_2)_{2:3,:}
    \end{bmatrix},
  \f]
  where \f$ (R_c R_i)_{2:3,:} \f$ denotes the last two rows of \f$ R_c R_i \f$. The
  rotation about the normal (spin) is not taken into account. The first component of
  \f$ y \f$ is the gap \f$ g \f$, to be set by the derived class (e.g.
  ContactVelocityAngularVelocityR in mechanics). The impulse applied to the bodies is
  \f$ H^\top \lambda \f$, for a NewtonImpactRollingFrictionNSL.

  Only the first row of \f$ J_h \f$ (the gradient of the gap, see NewtonEuler1DR) is
  computed; the other rows are zero.

  \see DevNotes, chapter "NewtonEulerR: computation of \f$ \nabla_q H \f$", subsection
  "Case FC3D: using the local frame local velocities".
*/
class NewtonEulerVelocityAngularVelocityR : public NewtonEuler1DR {
 private:
  ACCEPT_SERIALIZATION(NewtonEulerVelocityAngularVelocityR);

 public:
  /** destructor
   */
  virtual ~NewtonEulerVelocityAngularVelocityR() noexcept = default;

  /** allocate \f$ J_h \f$ (5 rows) and the operators of NewtonEulerR.
   *
   *  \param inter the interaction using this relation
   */
  void initialize(Interaction& inter) override;

  /** compute \f$ H(q) \f$ directly from the contact point and the contact frame (instead of
   *  the product \f$ J_h T \f$).
   *
   *  \param inter interaction that owns the relation
   *  \param q0  coordinates \f$ q \f$ of the bodies of the interaction
   */
  void computeH_NE_prod_T(const Interaction& inter,
                          const siconos::algebra::BlockVector& q0) override;

  virtual void accept(relations::Visitor& tourist) const override { tourist.visit(*this); }
};

}  // namespace siconos::modeling
#endif  // NEWTONEULERVELOCITYANGULARVELOCITYR_HPP
