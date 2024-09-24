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

/*! \file CablDS.hpp
  FEM model for Cable dynamical systems.

*/

#ifndef CABLEDS_H
#define CABLEDS_H

#include "LagrangianDS.hpp"

namespace siconos::fem::cable {
/**
   Cable-like dynamical systems

   Todo : document this class

  CableDS class - Second Order Non Linear Dynamical Systems, cable model after finite-element
  discretization.

  The equation of motion is also shortly denoted as  \f$ M(q,z) \dot v = F(v, q, t) + p \f$

  where  \f$ F(v, q, t) \in R^{ndof} \f$ collects the total forces acting on the
  system, that is \f$ F(v, q, t, z) =  F_{ext}(t) -  F_{int}(v, q , t) \f$


  This vector is saved and may be accessed using forces() method.

  Todo: check and update doc and comments.

  Add ref to Charlelie's phd report.

*/
class CableDS : public siconos::modeling::LagrangianDS {
 protected:
  double _EA{1};
  double _l_e{1};

  std::shared_ptr<siconos::algebra::SiconosMatrix> TRNp_Np{nullptr};

  void matmult(const std::shared_ptr<siconos::algebra::SiconosVector> &V, size_t a_startIdx,
               std::shared_ptr<siconos::algebra::SiconosVector> &R);
  void matmult2(const std::shared_ptr<siconos::algebra::SiconosVector> &V,
                std::shared_ptr<siconos::algebra::SiconosMatrix> &R);

  CableDS() = delete;
  CableDS(const CableDS &) = delete;
  CableDS(CableDS &&) = delete;
  CableDS &operator=(const CableDS &) = delete;
  CableDS &operator=(CableDS &&) = delete;

 public:
  CableDS(Eigen::Ref<siconos::algebra::SiconosVector> q0,
          Eigen::Ref<siconos::algebra::SiconosVector> velocity0,
          Eigen::Ref<siconos::algebra::SiconosMatrix> mass, double a_EA,
          double a_elem_length);

  ~CableDS() noexcept = default;

  // Methods that must be overloaded

  // This function will be called by the integrator at each time
  // step to update  \f$ F(v, q, t, z) \f$
  // --> takes into account fInt and fext
  void computeForces(double time, std::shared_ptr<siconos::algebra::SiconosVector> q,
                     std::shared_ptr<siconos::algebra::SiconosVector> velocity) override;

  // \f$ \nabla_q F \f$
  void computeJacobianqForces(double time) override;

  // \f$ \nabla_v F \f$
  void computeJacobianvForces(double time) override;

  void tangentStiffnessMatrix(std::shared_ptr<siconos::algebra::SiconosVector> q);
  void dampingMatrix();
  // + some access op to be added later, if required

  std::shared_ptr<siconos::algebra::SiconosMatrix> TRNp_NpMatrix();
};
}  // namespace siconos::fem::cable

#endif
