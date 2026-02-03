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
/*! \file NewtonEulerR.hpp

 */
#ifndef NEWTONEULERRELATION_H
#define NEWTONEULERRELATION_H

#include <cassert>
#include <optional>

#include "Relation.hpp"

namespace siconos::modeling {

class DynamicalSystem;

/**
   NewtonEuler (Non Linear) Relation (generic interface), with

   \f[
        y &=& H_NE.q + e \\
        \dot y &=& H_{NE}.Tv \\

        R = (H_{NE}.T)\top \lambda
   \f]

    Up to now, We are only considering holonomic relations.

  The following operators can be set by a user-defined function:

   - \f$ \frac{\partial }{\partial t} H_{NE} \f$

  Storage:

    - for \f$ H_NE \f$:
          - shared view for NewtonEulerR when setConstantH is called (i.e. no internal storage)
          - internal for all derived classes (since H is explicitely computed by a dedicated
  function)
    - for \f$ \frac{\partial }{\partial t}H_{NE} \f$ only at first call of
  computeSecondOrderTimeDerivativeTerms or if setComputeH_NE_dotFunction is called.
    - for e: shared view (i.e. no internal storage) if setConstanteVector is called
    - for \f$ H_NE . T\f$: always.
    - for T: always (in T_buffer_)

  All allocations are done during initialize call, except for \f$ \frac{\partial }{\partial
  t}H_{NE} \f$.
 */

class NewtonEulerR : public Relation {
 public:
  // add deltaq ??? -- xhub 30/03/2014
  /** Names and positions of ds variables that might be used (read only) in compute input and
   * output.
   * Last param, size, gives the size of the BlockVector to be allocated.
   */
  enum class ds_var : std::size_t { q0, velocity, dotq, p0, p1, p2, size };

 protected:
  ACCEPT_SERIALIZATION(NewtonEulerR);

  /** View onto the H_NE matrix of the relation */
  std::shared_ptr<siconos::algebra::MapType> H_NE_view_{nullptr};

  /** internal (optional) storage used for \f$ H_{NE} \f$ */
  siconos::algebra::SiconosMatrix H_NE_internal_storage_{0, 0};

  /** True if \f$ H_{NE} \f$ is a constant matrix */
  bool hasConstantH_NE_{false};

  // Rq: no need for internal storage. jacobianhOver_q is set with setHMatrix and is only
  // a view to some external storage.

  /** \f$ \frac{\partial}{\partial t}H_{NE} \f$
   *  This value is useful to compute the second-order
   *  derivative of the constraints with respect to time.
   */
  std::shared_ptr<siconos::algebra::SiconosMatrix> H_NE_dot_{nullptr};

  /** function wrapper used to compute \f$ \frac{\partial}{\partial t}H_{NE} \f$
   */
  siconos::modeling::func_prototypes::FunctionBVBV_M computeH_NE_dot_{nullptr};

  /** True if  \f$ \frac{\partial}{\partial t}H_{NE} \f$ is taken into account */
  bool hasH_NE_dot_{false};

  /** e */  // Used for the projection formulation
  std::shared_ptr<siconos::algebra::MapVectorType> eVector_view_{nullptr};

  /** True if e(t) is taken into account */
  bool haseVector_{false};

  /** vector of contact forces, ie: contactForce_ = B lambda. Useful for the end user.*/
  std::shared_ptr<siconos::algebra::SiconosVector> contactForce_{nullptr};

  /** Internal storage/buffer used to save H.T product.
     In the case of the bilateral constrains, it is H.T.
     In the case of a local frame, H.T is built from the geometrical
     datas(local frame, point of contact).*/
  std::shared_ptr<siconos::algebra::SiconosMatrix> H_NE_prod_T_{nullptr};

  /** buffer to save \f[ T(q) = \left[\begin{array}{cc} I_{3x3}  & 0 \\
               0 &  \phi(p) \end{array}\right] \f]
      from dynamical systems
  */
  std::shared_ptr<siconos::algebra::SiconosMatrix> T_buffer_{nullptr};

  /**  the additional  terms of the second order time derivative of y
   *
   *   \f$ \nabla_q h(q) \dot T v + \frac{d}{dt}H_{NE}  T v \f$
   *
   */
  std::shared_ptr<siconos::algebra::SiconosVector> secondOrderTimeDerivativeTerms_{nullptr};

  /** update H_NE matrix
   *
   *  \param time current time
   *  \param inter the interaction using this relation
   *  \param q0  q states vectors of the related the dynamical systems
   */
  virtual void computeH_NE_(double time, Interaction& inter,
                            const siconos::algebra::BlockVector& q0) {}

 public:
  /** Default and only constructor */
  NewtonEulerR() : Relation(RelationType::NewtonEuler, RelationSubType::NonLinearR) {};

  /** destructor */
  virtual ~NewtonEulerR() noexcept = default;

  /** initialize components specific to derived classes.
   *
   *  \param inter  Interaction associated with the Relation
   */
  virtual void initialize(Interaction& inter) override;

  /** check sizes of the relation specific operators.
   *
   *  \param inter an Interaction using this relation
   */
  virtual void checkSize(const Interaction& inter) const override;

  /** \return a read-only view on \f$ H(q, \ldots) \f$ matrix */
  inline auto H_NE() const {
    return siconos::algebra::ConstMapType(H_NE_view_->data(), H_NE_view_->rows(),
                                          H_NE_view_->cols());
  }

  /** Set a constant \f$ H(q, \ldots) \f$ matrix for the system (warning: shared memory)
   *
   *  \param newValue H matrix
   *
   */
  void setConstantH_NE(Eigen::Ref<siconos::algebra::SiconosMatrix> newValue);

  /** \return  a read-only view on e(t) */
  inline auto eVector() const {
    return siconos::algebra::ConstMapVectorType(eVector_view_->data(), eVector_view_->size());
  }

  /** set a constant e vector
   *
   *  \param neweVector e vector
   */
  void setConstanteVector(Eigen::Ref<siconos::algebra::SiconosVector> neweVector);

  /** True if e(t) is taken into account */
  bool haseVector() const { return haseVector_; }

  /** \return a read-only view on \f$ H(q, \ldots).T \f$ matrix */
  inline auto H_NE_prod_T() const {
    return siconos::algebra::ConstMapType(H_NE_prod_T_->data(), H_NE_prod_T_->rows(),
                                          H_NE_prod_T_->cols());
  }

  /** set a user-defined function to compute
   *  \f$ \frac{\partial}{\partial t}(\nabla^T_{q} h(q))\f$
   *
   *  \param fct the user-defined function (std::function, lambda ...)
   */
  void setComputeH_NE_dotFunction(
      const siconos::modeling::func_prototypes::FunctionBVBV_M& fct);

  /** Update \f$ \frac{\partial}{\partial t}(\nabla^T_{q} h(q))\f$
   *  \param q 'list' of state vectors (for all ds involved in the interaction)
   *  \param qdot 'list' of state vectors (for all ds involved in the interaction)
   */
  virtual void computeH_NE_dot(const siconos::algebra::BlockVector& q,
                               const siconos::algebra::BlockVector& qdot);

  /** \return  a read-only view on additional  terms of the second order time derivative of y
   *  \f$ \nabla_q h(q) \dot T v + \frac{d}{dt}(\nabla_q h(q) ) T v \f$
   */
  inline auto secondOrderTimeDerivativeTerms() {
    return siconos::algebra::ConstMapVectorType(secondOrderTimeDerivativeTerms_->data(),
                                                secondOrderTimeDerivativeTerms_->size());
  }

  /**
       To compute the terms of the second order time derivative of y
       \f$  \nabla_q h(q) \dot T v + \frac{d}{dt}(\nabla_q h(q) ) T v \f$

       \param time  current time
       \param inter interaction that owns the relation
       \param ds1 dynamical system linked to this interaction (source)
       \param ds2 second ds linked to this interaction (target). If there is
       only one ds in the inter, call this function with ..., ds, ds)
    */
  void computeSecondOrderTimeDerivativeTerms(double time, Interaction& inter,
                                             std::shared_ptr<DynamicalSystem> ds1,
                                             std::shared_ptr<DynamicalSystem> ds2);

  /** \return  a read-only view on the vector of the forces due to this relation
   */
  inline auto contactForce() const {
    return siconos::algebra::ConstMapVectorType(contactForce_->data(), contactForce_->size());
  }

  /**
      to compute the output y = h(q) of the Relation

      \param[in] q generalized coordinates vector of the concerned dynamical systems
      \param[in,out] y the resulting vector
  */
  virtual void computeh(const siconos::algebra::BlockVector& q,
                        Eigen::Ref<siconos::algebra::SiconosVector> y);

  /**
     to compute the output y = h(q) of the Relation

      \param[in] q1 generalized coordinates vector of the fist dynamical system involved
      in the relation
      \param[in] q2 generalized coordinates vector of the second dynamical system
      involved in the relation
      \param[in,out] y the resulting vector
  */
  virtual void computeh(
      const Eigen::Ref<const siconos::algebra::SiconosVector>& q1,
      const std::optional<Eigen::Ref<const siconos::algebra::SiconosVector>>& q2,
      Eigen::Ref<siconos::algebra::SiconosVector> y);

  /** default implementation consists in multiplying jachq and T
   *  in this implementation we use T_buffer_ which is consitent which directly
   *  computed with computeT(q) when q is given
   *  this one in more consistent with the notion of function of q
   *
   *  \param inter interaction that owns the relation
   *  \param q0  the block vector to the dynamical system position
   */
  virtual void computeH_NE_prod_T(const Interaction& inter,
                                  const siconos::algebra::BlockVector& q0);

  /** Update H_NE and \f$ H_NE. T \f$
   *
   *  \param time current time
   *  \param inter the interaction using this relation
   */
  virtual void computeJach(double time, Interaction& inter) override;

  /** to compute output
   *
   *  \param time current time
   *  \param inter the interaction using this relation
   *  \param derivativeNumber number of the derivative to compute, optional,
   *  default = 0.
   */
  virtual void computeOutput(double time, Interaction& inter,
                             unsigned int derivativeNumber = 0) override;

  /** to compute the input
   *
   *  \param time current time
   *  \param inter the interaction using this relation
   *  \param level number of the derivative to compute, optional, default = 0.
   */
  virtual void computeInput(double time, Interaction& inter, unsigned int level = 0) override;

  void display() const override {}

  /** @brief allocation of memory space for relations in the graph
   *      Warning: internal use only (called from Topology)
   *  @param dslink a container of vectors (pointers), from the parent interaction
   *  @param ds1 first ds concerned by the relation
   *  @param ds2 second ds concerned by the relation
   */
  virtual void allocate_read_dynamical_systems_var_vectors(
      std::vector<std::shared_ptr<siconos::algebra::BlockVector>>& ds_vars,
      modeling::DynamicalSystem& ds1, modeling::DynamicalSystem& ds2) const override;

  virtual void accept(relations::Visitor& tourist) const override { tourist.visit(*this); }
};
}  // namespace siconos::modeling
#endif  // NEWTONEULERRELATION_H
