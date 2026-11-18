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
/*! \file NewtonEuler3DR.hpp
 *  \brief Contact relation between Newton-Euler bodies: relative velocity in the contact
 *  frame (normal and tangential components, for Coulomb friction).
 */
#ifndef NEWTONEULER3DR_HPP
#define NEWTONEULER3DR_HPP

#include "NewtonEuler1DR.hpp"
#include "SiconosMatrix.hpp"

namespace siconos::modeling {

/** Unilateral contact with friction between body 1 and body 2 (or body 1 and a fixed
  obstacle): relative velocity at the contact point, in the contact frame.

  The geometry (contact points \f$ P_1, P_2 \f$, normal \f$ n \f$, gap \f$ g \f$), the
  contact frame \f$ R_c = [\, n \;\; t_1 \;\; t_2 \,]^\top \f$ and the matrix \f$ H_c(q) \f$
  are defined in NewtonEuler1DR. Here \f$ m = 3 \f$ and
  \f[
    \dot y = \begin{bmatrix} u_N \\ u_{T_1} \\ u_{T_2} \end{bmatrix} = R_c\, u = H_c(q)\, v,
    \qquad
    H_c(q) = R_c \begin{bmatrix}
        I_3 & [x_{G_1} - P]_\times R_1 & -I_3 & -[x_{G_2} - P]_\times R_2
      \end{bmatrix},
  \f]
  with \f$ u \f$ the relative velocity of body 1 with respect to body 2 at
  \f$ P = P_1 \f$. The first component of \f$ y \f$ is the gap \f$ g \f$, to be set by the
  derived class (e.g. ContactR in mechanics). The impulse applied to the bodies is
  \f$ H_c^\top \lambda \f$, with \f$ \lambda = (r_N, r_{T_1}, r_{T_2}) \f$ the contact
  impulse in the contact frame, for a NewtonImpactFrictionNSL.

  Only the first row of \f$ J_h \f$ (the gradient of the gap, see NewtonEuler1DR) is
  computed; the other rows are zero.

  \see DevNotes, chapter "NewtonEulerR: computation of \f$ \nabla_q H \f$", subsection
  "Case FC3D: using the local frame local velocities".
*/
class NewtonEuler3DR : public NewtonEuler1DR {
 private:
  ACCEPT_SERIALIZATION(NewtonEuler3DR);

 public:
  /** destructor */
  virtual ~NewtonEuler3DR() noexcept = default;

  /** allocate \f$ J_h \f$ (3 rows) and the operators of NewtonEulerR.
   *
   *  \param inter the interaction using this relation
   */
  void initialize(Interaction& inter) override;

  /** compute \f$ H = H_c(q) \f$ directly from the contact point and the contact frame
   *  (instead of the product \f$ J_h T \f$).
   *
   *  \param inter  interaction that owns the relation
   *  \param q0  coordinates \f$ q \f$ of the bodies of the interaction
   */
  void computeH_NE_prod_T(const Interaction& inter,
                          const siconos::algebra::BlockVector& q0) override;

  void display() const override {}
  virtual void accept(relations::Visitor& tourist) const override { tourist.visit(*this); }
};
}  // namespace siconos::modeling
#endif  // NEWTONEULER3DR_HPP
