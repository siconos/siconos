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
/*! \file Lagrangian2d#DR.hpp

 */
#ifndef NodeFem2d3DR_H
#define NodeFem2d3DR_H
#include "LagrangianScleronomousR.hpp"

/** NodeFem2d3DR
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
class NodeFem2d3DR : public siconos::modeling::LagrangianScleronomousR {
 protected:
  //  ACCEPT_SERIALIZATION(NodeFem2d3DR);

  /* index of the node of the Fem cable involved in this relation */
  unsigned int _node_index;

  /* Current Contact Points */
  siconos::algebra::SiconosVector3 contactPoint1_ = siconos::algebra::SiconosVector3::Zero();
  siconos::algebra::SiconosVector3 contactPoint2_ = siconos::algebra::SiconosVector3::Zero();

  /* Normal vector at contact */
  siconos::algebra::SiconosVector3 nc_ = siconos::algebra::SiconosVector3::Zero();

  /** Tangent vector at contact */
  siconos::algebra::SiconosVector3 tangent_ = siconos::algebra::SiconosVector3::Zero();

  /** Set the coordinates of first contact point.  Must only be done
   *  in a computeh() override.
   *
   *  \param npc new coordinates
   */
  void setpc1(const siconos::algebra::SiconosVector3 &npc) { contactPoint1_ = npc; };

  /** Set the coordinates of second contact point.  Must only be done
   *  in a computeh() override.
   *
   *  \param npc new coordinates
   */
  void setpc2(const siconos::algebra::SiconosVector3 &npc) { contactPoint2_ = npc; };

  /** Set the coordinates of inside normal vector at the contact point.

   *  Must only be done in a computeh() override.
   *
   *  \param nnc new coordinates
   */
  void setnc(const siconos::algebra::SiconosVector3 &nnc) { nc_ = nnc; };

  /** Set the coordinates of inside normal vector at the contact point.
   *  Must only be done in a computeh() override.
   *
   *  \param nnc new coordinates
   */
  void settc(const siconos::algebra::SiconosVector3 &ntc) { tangent_ = ntc; };

 public:
  /** constructor
   */
  NodeFem2d3DR(unsigned int node_index) : LagrangianScleronomousR(), _node_index(node_index) {}

  /** constructor
   */
  NodeFem2d3DR(unsigned int node_index, const siconos::algebra::SiconosVector3 &pc1,
               const siconos::algebra::SiconosVector3 &pc2,
               const siconos::algebra::SiconosVector3 &normal,
               const siconos::algebra::SiconosVector3 &tangent)
      : LagrangianScleronomousR(),
        _node_index(node_index),
        contactPoint1_(pc1),
        contactPoint2_(pc2),
        nc_(normal),
        tangent_(tangent) {}

  /** destructor
   */
  virtual ~NodeFem2d3DR() noexcept {};

  void initialize(siconos::modeling::Interaction &inter) override;

  /**
     to compute the output y = h(q,z) of the Relation

     \param q coordinates of the dynamical systems involved in the relation
     \param z user defined parameters (optional)
     \param y the resulting vector
  */
  void computeh(const siconos::algebra::BlockVector &q,
                Eigen::Ref<siconos::algebra::SiconosVector> y) override;
  /**
     to compute the jacobian of h(...). Set attribute _jachq (access: jacqhq())

     \param q coordinates of the dynamical systems involved in the relation
     \param z user defined parameters (optional)
  */
  void computeJacobianhOver_q(const siconos::algebra::BlockVector &q) override;

  /** Return the distance between pc1 and pc, with sign according to normal */
  double distance() const;

  inline const siconos::algebra::SiconosVector3 &pc1() const { return contactPoint1_; }

  inline const siconos::algebra::SiconosVector3 &pc2() const { return contactPoint2_; }

  inline const siconos::algebra::SiconosVector3 &nc() const { return nc_; }
  inline const siconos::algebra::SiconosVector3 &tangent() const { return tangent_; }

  /** update the contact points from references
   */
  void updateContactPoints(const siconos::algebra::SiconosVector3 &pc1,
                          const siconos::algebra::SiconosVector3 &pc2,
                          const siconos::algebra::SiconosVector3 &normal,
                          const siconos::algebra::SiconosVector3 &tangent) {
    contactPoint1_ = pc1;
    contactPoint2_ = pc2;
    nc_ = normal;
    tangent_ = tangent;
  };

  /** update the contact points from array
   */
  void updateContactPoints(double pc1[3], double pc2[3], double normal[3], double tangent[3]) {
    contactPoint1_ << pc1[0], pc1[1], pc1[2];
    contactPoint2_ << pc2[0], pc2[1], pc2[2];
    nc_ << normal[0], normal[1], normal[2];
    tangent_ << tangent[0], tangent[1], tangent[2];
  };
  void display() const override;

  // ACCEPT_STD_VISITORS();
};
}  // namespace siconos::mechanics::fem
#endif  // NodeFem2d3DR_H
