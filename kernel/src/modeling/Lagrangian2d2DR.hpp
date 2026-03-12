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
/*! \file Lagrangian2d2DR.hpp

 */
#ifndef Lagrangian2d2DR_H
#define Lagrangian2d2DR_H

#include "LagrangianScleronomousR.hpp"

namespace siconos::modeling {
/** Lagrangian2d2DR
 *
 * This class is an interface for a relation with impact.  It
 * implements the computation of the jacoboian of h from the points of
 * contacts and the normal.  Use this class consists in overloading
 * the method computeh, by setting the member pc1, pc2, nc and y.  The
 * matrix jachq is used both for the building of the OSNSP (with T)
 * and for the predictor of activation of deactivation of the Interaction.
 *
 */

class Lagrangian2d2DR : public LagrangianScleronomousR {
 protected:
  ACCEPT_SERIALIZATION(Lagrangian2d2DR);

  /* Current Contact Points, may be updated within Newton loop based
   * on _relPc1, _relPc2. */
  std::shared_ptr<siconos::algebra::SiconosVector> _Pc1{nullptr};
  std::shared_ptr<siconos::algebra::SiconosVector> _Pc2{nullptr};

  /* Inward Normal at the contact.
   * \todo The meaning of "Inward" has to be explained carefully.
   */
  std::shared_ptr<siconos::algebra::SiconosVector> _Nc{nullptr};

  /** parametrize the projection algorithm
   *  the projection is done on the surface \f$y=0\f$ or on \f$y \geq 0\f$
   */
  bool _isOnContact = false;

 public:
  /** constructor */
  Lagrangian2d2DR()
      : LagrangianScleronomousR(),
        _Pc1{std::make_shared<siconos::algebra::SiconosVector>(2)},
        _Pc2{std::make_shared<siconos::algebra::SiconosVector>(2)},
        _Nc{std::make_shared<siconos::algebra::SiconosVector>(2)} {
    /*_ds1=nullptr;_ds2=nullptr;*/
  }

  /** destructor
   */
  virtual ~Lagrangian2d2DR() noexcept = default;

  void initialize(Interaction &inter) override;

  /** to compute the output y = h(q,z) of the Relation
     \param q coordinates of the dynamical systems involved in the relation
     \param z user defined parameters (optional)
     \param y the resulting vector
 */
  void computeh(const siconos::algebra::BlockVector &q,
                Eigen::Ref<siconos::algebra::SiconosVector> y) override;

  /** to compute the jacobian of h(...). Set attribute jacobianhOver_q_ (access: jacqhq())
      \param q coordinates of the dynamical systems involved in the relation
      \param z user defined parameters (optional)
  */
  void computeJacobianhOver_q(const siconos::algebra::BlockVector &q) override;

  /** Return the distance between pc1 and pc, with sign according to normal */
  double distance() const;

  inline auto pc1() const { return _Pc1; }
  inline auto pc2() const { return _Pc2; }
  inline auto nc() const { return _Nc; }
  void display() const override;
  virtual void accept(relations::Visitor &tourist) const override { tourist.visit(*this); }
};
}  // namespace siconos::modeling
#endif  // NEWTONEULERRIMPACT_H
