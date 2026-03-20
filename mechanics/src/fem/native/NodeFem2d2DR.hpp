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
/*! \file Lagrangian2d#DR.hpp

 */
#ifndef NodeFem2d2DR_H
#define NodeFem2d2DR_H

#include "LagrangianScleronomousR.hpp"

/** NodeFem2d2DR
 *
 * This class is an interface for a relation with impact.  It
 * implements the computation of the jacoboian of h from the points of
 * contacts and the normal.  Use this class consists in overloading
 * the method computeh, by setting the member pc1, pc2, nc and y.  The
 * matrix jachq is used both for the building of the OSNSP
 * and for the predictor of activation of deactivation of the Interaction.
 *
 */
namespace siconos::mechanics::fem {

class FENode;

class NodeFem2d2DR : public modeling::LagrangianScleronomousR {
 protected:
  /** index of the node of the Fem cable involved in this relation */
  std::shared_ptr<FENode> node_{nullptr};

  /* Current Contact Points */
  siconos::algebra::SiconosVector2 contactPoint1_ = siconos::algebra::SiconosVector2::Zero();
  siconos::algebra::SiconosVector2 contactPoint2_ = siconos::algebra::SiconosVector2::Zero();

  /* Normal vector at contact */
  siconos::algebra::SiconosVector2 nc_ = siconos::algebra::SiconosVector2::Zero();

  /** Tangent vector at contact */
  siconos::algebra::SiconosVector2 tangent_ = siconos::algebra::SiconosVector2::Zero();

 public:
  /** constructor
   */
  NodeFem2d2DR(std::shared_ptr<FENode> node) : LagrangianScleronomousR(), node_(node) {}

  /** constructor
   */
  NodeFem2d2DR(std::shared_ptr<FENode> node, const siconos::algebra::SiconosVector2 &pc2,
               const siconos::algebra::SiconosVector2 &normal,
               const siconos::algebra::SiconosVector2 &tangent)
      : LagrangianScleronomousR(),
        node_(node),
        contactPoint2_(pc2),
        nc_(normal),
        tangent_(tangent) {}

  /** destructor
   */
  virtual ~NodeFem2d2DR() noexcept = default;

  void initialize(siconos::modeling::Interaction &inter) override;

  /**
     to compute the output y = h(q,z) of the Relation

     \param q coordinates of the dynamical systems involved in the relation
     \param z user defined parameters (optional)
     \param y the resulting vector
  */
  void computeh(const siconos::algebra::BlockVector &q,
                Eigen::Ref<siconos::algebra::SiconosVector> y) override;

  /** Return the distance between pc1 and pc, with sign according to normal */
  double distance() const;

  inline std::shared_ptr<FENode> node() const { return node_; }
  inline const siconos::algebra::SiconosVector2 &pc1() const { return contactPoint1_; }

  inline const siconos::algebra::SiconosVector2 &pc2() const { return contactPoint2_; }

  inline const siconos::algebra::SiconosVector2 &nc() const { return nc_; }
  inline const siconos::algebra::SiconosVector2 &tangent() const { return tangent_; }

  /** update the contact points from references
   */
  void updateContactPoint(const siconos::algebra::SiconosVector2 &pc2,
                          const siconos::algebra::SiconosVector2 &normal,
                          const siconos::algebra::SiconosVector2 &tangent);

  /** update the contact points from array
   */
  void updateContactPoint(double pc2[2], double normal[2], double tangent[2]);

  void display() const override;

  // ACCEPT_STD_VISITORS();
};
}  // namespace siconos::mechanics::fem
#endif  // NodeFem2d2DR_H
