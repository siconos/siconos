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
  std::shared_ptr<siconos::algebra::SiconosVector> _Pc1{nullptr};
  std::shared_ptr<siconos::algebra::SiconosVector> _Pc2{nullptr};

  /* Normal vector at contact.
   */
  std::shared_ptr<siconos::algebra::SiconosVector> _Normal{nullptr};

  /* Tangent vector at contact.
   */
  std::shared_ptr<siconos::algebra::SiconosVector> _Tangent{nullptr};

  /** Set the coordinates of first contact point.  Must only be done
   *  in a computeh() override.
   *
   *  \param npc new coordinates
   */
  void setpc1(std::shared_ptr<siconos::algebra::SiconosVector> npc) { _Pc1 = npc; };

  /** Set the coordinates of second contact point.  Must only be done
   *  in a computeh() override.
   *
   *  \param npc new coordinates
   */
  void setpc2(std::shared_ptr<siconos::algebra::SiconosVector> npc) { _Pc2 = npc; };

  /** Set the coordinates of inside normal vector at the contact point.

   *  Must only be done in a computeh() override.
   *
   *  \param nnc new coordinates
   */
  void setnc(std::shared_ptr<siconos::algebra::SiconosVector> nnc) { _Normal = nnc; };

  /** Set the coordinates of inside normal vector at the contact point.
   *  Must only be done in a computeh() override.
   *
   *  \param nnc new coordinates
   */
  void settc(std::shared_ptr<siconos::algebra::SiconosVector> ntc) { _Tangent = ntc; };

 public:
  /** constructor
   */
  Cable2d3DR(unsigned int node_index)
      : LagrangianScleronomousR(),
        _node_dof_index(node_index),
        _Pc1{std::make_shared<siconos::algebra::SiconosVector>(3)},
        _Pc2{std::make_shared<siconos::algebra::SiconosVector>(3)},
        _Normal{std::make_shared<siconos::algebra::SiconosVector>(3)},
        _Tangent{std::make_shared<siconos::algebra::SiconosVector>(3)} {}

  /** constructor
   */
  Cable2d3DR(unsigned int node_index, std::shared_ptr<siconos::algebra::SiconosVector> pc2,
             std::shared_ptr<siconos::algebra::SiconosVector> normal,
             std::shared_ptr<siconos::algebra::SiconosVector> tangent)
      : LagrangianScleronomousR(),
        _node_dof_index(node_index),
        _Pc1{std::make_shared<siconos::algebra::SiconosVector>(3)},
        _Pc2(pc2),
        _Normal(normal),
        _Tangent(tangent) {}

  /** destructor */
  virtual ~Cable2d3DR() noexcept = default;

  void initialize(siconos::modeling::Interaction &inter) override;

  /**
     to compute the output y = h(q,z) of the Relation

     \param q coordinates of the dynamical systems involved in the relation
     \param z user defined parameters (optional)
     \param y the resulting vector
  */
  void computeh(const siconos::algebra::BlockVector &q, siconos::algebra::BlockVector &z,
                siconos::algebra::SiconosVector &y) override;

  /**
     to compute the jacobian of h(...). Set attribute _jachq (access: jacqhq())

     \param q coordinates of the dynamical systems involved in the relation
     \param z user defined parameters (optional)
  */
  void computeJachq(const siconos::algebra::BlockVector &q,
                    siconos::algebra::BlockVector &z) override;

  /** Return the distance between pc1 and pc, with sign according to normal */
  double distance() const;

  inline std::shared_ptr<siconos::algebra::SiconosVector> pc1() const { return _Pc1; }
  inline std::shared_ptr<siconos::algebra::SiconosVector> pc2() const { return _Pc2; }
  inline std::shared_ptr<siconos::algebra::SiconosVector> normal() const { return _Normal; }
  inline std::shared_ptr<siconos::algebra::SiconosVector> tangent() const { return _Tangent; }

  /** update the contact points
   * this method is a bit useless if the relation has been constructed with
   * shared pointer
   */
  void updateContactPoint(std::shared_ptr<siconos::algebra::SiconosVector> pc1,
                          std::shared_ptr<siconos::algebra::SiconosVector> pc2,
                          std::shared_ptr<siconos::algebra::SiconosVector> normal,
                          std::shared_ptr<siconos::algebra::SiconosVector> tangent);

  /** update the contact points from references
   */
  void updateContactPoint(siconos::algebra::SiconosVector &pc1,
                          siconos::algebra::SiconosVector &pc2,
                          siconos::algebra::SiconosVector &normal,
                          siconos::algebra::SiconosVector &tangent);

  /** update the contact points from array
   */
  void updateContactPoint(double pc1[3], double pc2[3], double normal[3], double tangent[3]);

  void display() const override;
};
}  // namespace siconos::fem::cable
#endif  // Cable2d3DR_H
