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
/*! \file Cable2d3DR.h


 */
#ifndef Cable2d3DR_H
#define Cable2d3DR_H

#include "LagrangianScleronomousR.hpp"

namespace siconos::fem::cable {

/** Cable2d3DR
 *
 * This class is an interface for a relation with impact.  It
 * implements the computation of the jacoboian of h from the points of
 * contacts and the normal.  Use this class consists in overloading
 * the method computeh, by setting the member pc1, pc2, nc and y.  The
 * matrix jachq is used both for the building of the OSNSP
 * and for the predictor of activation of deactivation of the Interaction.
 *
 */
class Cable2d3DR : public siconos::modeling::LagrangianScleronomousR {
 protected:
  /* index of the node of the Fem cable involved in this relation */
  unsigned int _node_dof_index{0};

  /* Current Contact Points, may be updated within Newton loop based
   * on _relPc1, _relPc2. */
  siconos::algebra::SiconosVector3 contactPoint1_ = siconos::algebra::SiconosVector3::Zero();
  siconos::algebra::SiconosVector3 contactPoint2_ = siconos::algebra::SiconosVector3::Zero();

  /**  Normal vector at contact */
  siconos::algebra::SiconosVector3 nc_ = siconos::algebra::SiconosVector3::Zero();

  /** Tangent vector at contact */
  siconos::algebra::SiconosVector3 tangent_ = siconos::algebra::SiconosVector3::Zero();

 public:
  /** constructor */
  Cable2d3DR(unsigned int node_index)
      : LagrangianScleronomousR(), _node_dof_index{node_index} {}

  /** constructor */
  Cable2d3DR(unsigned int node_index, const siconos::algebra::SiconosVector3 &pc2,
             const siconos::algebra::SiconosVector3 &normal,
             const siconos::algebra::SiconosVector3 &tangent)
      : LagrangianScleronomousR(),
        _node_dof_index{node_index},
        contactPoint2_(pc2),
        nc_(normal),
        tangent_(tangent) {}

  /** destructor */
  virtual ~Cable2d3DR() noexcept = default;

  void initialize(siconos::modeling::Interaction &inter) override;

  /**
   to compute the output y = h(q) of the Relation

   \param q coordinates of the dynamical systems involved in the relation
   \param y the resulting vector
 */
  void computeh(const siconos::algebra::BlockVector &q,
                Eigen::Ref<siconos::algebra::SiconosVector> y) override;

  /** Computes \f$ \nabla^\top_q h(q) \f$
   * \param q coordinates of the dynamical systems involved in the relation
   */
  void computeJacobianhOver_q(const siconos::algebra::BlockVector &q) override;

  /** Return the distance between pc1 and pc, with sign according to normal */
  double distance() const;
  inline const siconos::algebra::SiconosVector3 &pc1() const { return contactPoint1_; }
  inline const siconos::algebra::SiconosVector3 &pc2() const { return contactPoint2_; }
  inline const siconos::algebra::SiconosVector3 &normal() const { return nc_; }
  inline const siconos::algebra::SiconosVector3 &tangent() const { return tangent_; }

  /** update the contact points from references
   */
  void updateContactPoints(const siconos::algebra::SiconosVector3 &pc1,
                          const siconos::algebra::SiconosVector3 &pc2,
                          const siconos::algebra::SiconosVector3 &normal,
                          const siconos::algebra::SiconosVector3 &tangent);

  /** update the contact points from array
   */
  void updateContactPoints(double pc1[3], double pc2[3], double normal[3], double tangent[3]);

  void display() const override;
};
}  // namespace siconos::fem::cable
#endif  // Cable2d3DR_H
