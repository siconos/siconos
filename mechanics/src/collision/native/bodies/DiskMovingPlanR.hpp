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

/*! \file DiskMovingPlanR.hpp
 */
/**
   disk - moving plan relation - Inherits from LagrangianRheonomousR
*/

#ifndef DiskMovingPlanR_h
#define DiskMovingPlanR_h

#include <memory>

#include "LagrangianRheonomousR.hpp"
#include "PluginTypes.hpp"  // for FTime

namespace siconos::collision::native::bodies {
class DiskMovingPlanR : public siconos::modeling::LagrangianRheonomousR,
                        public std::enable_shared_from_this<DiskMovingPlanR> {
 private:
  ACCEPT_SERIALIZATION(DiskMovingPlanR);

  double _time{0.}, _A{0.}, _B{0.}, _C{0.}, _ADot{0.}, _BDot{0.}, _CDot{0.}, _sqrA2pB2{0.},
      _r{0.}, _AADot{0.}, _BBDot{0.}, _cubsqrA2pB2{0.};

  std::shared_ptr<siconos::plugins::PluggedObject> _AFunction{nullptr};
  std::shared_ptr<siconos::plugins::PluggedObject> _BFunction{nullptr};
  std::shared_ptr<siconos::plugins::PluggedObject> _CFunction{nullptr};

  std::shared_ptr<siconos::plugins::PluggedObject> _ADotFunction{nullptr};
  std::shared_ptr<siconos::plugins::PluggedObject> _BDotFunction{nullptr};
  std::shared_ptr<siconos::plugins::PluggedObject> _CDotFunction{nullptr};

  DiskMovingPlanR() = delete;

 public:
  DiskMovingPlanR(siconos::plugins::FTime, siconos::plugins::FTime, siconos::plugins::FTime,
                  siconos::plugins::FTime, siconos::plugins::FTime, siconos::plugins::FTime,
                  double);

  ~DiskMovingPlanR() noexcept = default;

  void init(double);

  /**
     to compute the output y = h(t,q,z) of the Relation

     \param time current time value
     \param q coordinates of the dynamical systems involved in the relation
     \param z user defined parameters (optional)
     \param y the resulting vector
  */
  void computeh(double time, const siconos::algebra::BlockVector& q,
                siconos::algebra::BlockVector& z, siconos::algebra::SiconosVector& y) override;

  /**
     to compute the jacobian of h(...). Set attribute _jachq (access: jacqhq())

     \param time current time value
     \param q coordinates of the dynamical systems involved in the relation
     \param z user defined parameters (optional)
  */
  void computeJachq(double time, const siconos::algebra::BlockVector& q,
                    siconos::algebra::BlockVector& z) override;

  /**
     to compute the time-derivative of the output y = h(t,q,z), saved in attribute _hDot
     (access: hDot())

     \param time current time value
     \param q coordinates of the dynamical systems involved in the relation
     \param z user defined parameters (optional)
  */
  void computehDot(double time, const siconos::algebra::BlockVector& q,
                   siconos::algebra::BlockVector& z) override;

  double distance(double, double, double);

  void setComputeAFunction(siconos::plugins::FTime f);

  void setComputeBFunction(siconos::plugins::FTime f);

  void setComputeCFunction(siconos::plugins::FTime f);

  void setComputeADotFunction(siconos::plugins::FTime f);

  void setComputeBDotFunction(siconos::plugins::FTime f);

  void setComputeCDotFunction(siconos::plugins::FTime f);

  bool equal(siconos::plugins::FTime, siconos::plugins::FTime, siconos::plugins::FTime,
             double) const;

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
};
}  // namespace siconos::collision::native::bodies

#endif /* DiskMovingPlanR */
