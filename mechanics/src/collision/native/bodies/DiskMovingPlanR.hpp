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

/*! \file DiskMovingPlanR.hpp
 */
/**
   disk - moving plan relation - Inherits from LagrangianRheonomousR
*/

#ifndef DiskMovingPlanR_h
#define DiskMovingPlanR_h

#include <memory>

#include "LagrangianRheonomousR.hpp"

namespace siconos::collision::native::bodies {
class DiskMovingPlanR : public siconos::modeling::LagrangianRheonomousR,
                        public std::enable_shared_from_this<DiskMovingPlanR> {
 private:
  ACCEPT_SERIALIZATION(DiskMovingPlanR);

  double _time{0.}, A_{0.}, B_{0.}, C_{0.}, Adot_{0.}, Bdot_{0.}, Cdot_{0.}, sqrA2pB2_{0.},
      radius_{0.}, AAdot_{0.}, BBdot_{0.}, cubsqrA2pB2_{0.};

  siconos::modeling::func_prototypes::FunctionS_S computeA_{nullptr};
  siconos::modeling::func_prototypes::FunctionS_S computeB_{nullptr};
  siconos::modeling::func_prototypes::FunctionS_S computeC_{nullptr};
  siconos::modeling::func_prototypes::FunctionS_S computeAdot_{nullptr};
  siconos::modeling::func_prototypes::FunctionS_S computeBdot_{nullptr};
  siconos::modeling::func_prototypes::FunctionS_S computeCdot_{nullptr};

  DiskMovingPlanR() = delete;

 public:
  /** default and only constructor
   * \param rad radius
   */
  DiskMovingPlanR(double rad);

  ~DiskMovingPlanR() noexcept = default;

  void init(double);

  /**
    to compute the output y = h(q, t) of the Relation

    \param q coordinates of the dynamical systems involved in the relation
    \param time current time value
    \param y the resulting vector
  */
  void computeh(const siconos::algebra::BlockVector &q, double time,
                Eigen::Ref<siconos::algebra::SiconosVector> y) override;

  /** Computes \f$ \nabla^\top_q h(q, t) \f$
   *  \param q coordinates of the dynamical systems involved in the relation
   *  \param time current time value
   */
  void computeJacobianhOver_q(const siconos::algebra::BlockVector &q, double time) override;

  /** Update \f$ \frac{\partial }{\partial t}h(q,t) \f$
   *  \param position 'list' of state vectors (for all ds involved in the interaction)
   *  \param time the current time
   */
  void computehdot(const siconos::algebra::BlockVector &position, double time) override;

  double distance(double, double, double);

  /** set a user-defined function to compute A(t)
   *
   *  \param fct the user-defined function (std::function, lambda ...)
   */
  void setComputeAFunction(const siconos::modeling::func_prototypes::FunctionS_S &fct);

  /** set a user-defined function to compute B(t)
   *
   *  \param fct the user-defined function (std::function, lambda ...)
   */
  void setComputeBFunction(const siconos::modeling::func_prototypes::FunctionS_S &fct);

  /** set a user-defined function to compute C(t)
   *
   *  \param fct the user-defined function (std::function, lambda ...)
   */
  void setComputeCFunction(const siconos::modeling::func_prototypes::FunctionS_S &fct);

  /** set a user-defined function to compute \f$ \frac{\partial }{\partial t} A(t) \f$
   *
   *  \param fct the user-defined function (std::function, lambda ...)
   */
  void setComputeAdotFunction(const siconos::modeling::func_prototypes::FunctionS_S &fct);

  /** set a user-defined function to compute \f$ \frac{\partial }{\partial t} B(t) \f$
   *
   *  \param fct the user-defined function (std::function, lambda ...)
   */
  void setComputeBdotFunction(const siconos::modeling::func_prototypes::FunctionS_S &fct);

  /** set a user-defined function to compute \f$ \frac{\partial }{\partial t} C(t) \f$
   *
   *  \param fct the user-defined function (std::function, lambda ...)
   */
  void setComputeCdotFunction(const siconos::modeling::func_prototypes::FunctionS_S &fct);

  bool equal(const siconos::modeling::func_prototypes::FunctionS_S &pA,
             const siconos::modeling::func_prototypes::FunctionS_S &pB,
             const siconos::modeling::func_prototypes::FunctionS_S &pC, double) const;

  /**
     compute A

     \param t the time
  */
  void computeA(double t);

  /**
     compute B

     \param t the time
  */
  void computeB(double t);

  /**
     compute C

     \param t the time
  */
  void computeC(double t);

  /**
     compute ADot

     \param t the time
  */
  inline void computeADot(double t);

  /**
     compute BDot

     \param t the time
  */
  inline void computeBDot(double t);

  /**
     compute CDot

     \param t the time
  */
  inline void computeCDot(double t);

  virtual void accept(modeling::relations::Visitor &tourist) const override { tourist.visit(*this); }

};
}  // namespace siconos::collision::native::bodies

#endif /* DiskMovingPlanR */
