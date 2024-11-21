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

/*! \file FirstOrderType1R.hpp
  non linear relations, with y depending on dynamical systems state and r
  on lambda.
*/

#ifndef FirstOrderType1R_H
#define FirstOrderType1R_H

#include "FirstOrderR.hpp"

namespace siconos::modeling {

/**
    First order non linear relations, with:

    \f[

       y &= h(X)\\
       R &= g(\lambda)
    \f]

    where X and R corresponds to DynamicalSystem variables.
    If 2 DynamicalSystem are involved in the Interaction, then X = [x1
    x2], R = [r1 r2].

    The following operators can be can be set by user-defined functions:

    - h(x)
    - \f$ g(\lambda) \f$
    - \f$ \nabla_{\lambda} g(\lambda) \f$
    - \f$ \nabla_{x} h(x) \f$

    Storage for:

    - h: in y vector from the interaction
    - g: in r vectors from the dynamical systems
    - \f$ \nabla_{x} h(x) \f$ in FirstOrderR base class matrix attribute
    - \f$ \nabla_{\lambda} g(\lambda) \f$ in FirstOrderR base class matrix attribute

*/
class FirstOrderType1R : public FirstOrderR {
 protected:
  ACCEPT_SERIALIZATION(FirstOrderType1R);

  /** function wrapper used to compute \f$ h(x,t,\lambda) \f$ */
  siconos::modeling::func_prototypes::FunctionBV_V computeh_{nullptr};

  /** function wrapper used to compute \f$ h(x,t,\lambda) \f$ */
  siconos::modeling::func_prototypes::FunctionV_BV computeg_{nullptr};

  /** function wrapper used to compute  \f$ \nabla_x h(x,t,\lambda) \f$ */
  siconos::modeling::func_prototypes::FunctionBV_M computejacobianhOver_state_{nullptr};

  /** function wrapper used to compute  \f$ \nabla_{\lambda} g(\lambda) \f$ */
  siconos::modeling::func_prototypes::FunctionV_M computejacobiangOver_lambda_{nullptr};

 public:
  /** default constructor */
  FirstOrderType1R() : FirstOrderR(RelationSubType::Type1R) {};

  /** destructor */
  virtual ~FirstOrderType1R() noexcept = default;

  /** initialize the relation (check sizes, memory allocation ...)
   *
   *  \param inter the interaction that owns this relation
   */
  void initialize(Interaction &inter) override;

  /** check sizes of the relation specific operators.
   *
   *  \param inter an Interaction using this relation
   */
  inline void checkSize(Interaction &inter) override;

  /** set a user-defined function to compute \f$ h(x) \f$
   *
   *  \param fct the user-defined function (std::function, lambda ...)
   */
  void setComputehFunction(const siconos::modeling::func_prototypes::FunctionBV_V &fct);

  /** To compute  \f$ h(x) \f$
   *
   *  \param[in] x state vector (for all DS involved in the relation)
   *  \param[out] y result, \f$ y\f$ value from interaction
   */
  virtual void computeh(const siconos::algebra::BlockVector &state,
                        Eigen::Ref<siconos::algebra::SiconosVector> y);

  /** set a user-defined function to compute \f$ g(\lambda) \f$
   *
   *  \param fct the user-defined function (std::function, lambda ...)
   */
  void setComputegFunction(const siconos::modeling::func_prototypes::FunctionV_BV &fct);

  /** To compute  \f$ g(\lambda) \f$
   *
   *  \param[in] lambda \f$ \lambda\f$ value from interaction
   *  \param[out] res result, \f$ r\f$ value (for all DS involved in the relation)
   */
  virtual void computeg(const Eigen::Ref<const siconos::algebra::SiconosVector> &lambda,
                        siconos::algebra::BlockVector &res);

  /** Set a constant \f$ \nabla_x h(x) \f$ matrix for the system
   *
   *  \param newValue \f$ \nabla_x h(x) \f$ matrix
   *
   */
  void setConstantJacobianhOver_state(Eigen::Ref<siconos::algebra::SiconosMatrix> newValue);

  /** set a user-defined function to compute \f$ \nabla_x h(x) \f$ \f$
   *
   *  \param fct the user-defined function (std::function, lambda ...)
   */
  void setComputeJacobianhOver_stateFunction(
      const siconos::modeling::func_prototypes::FunctionBV_M &fct);

  /** Computes \f$ \nabla_x h(x) \f$
   *  \param q coordinates of the dynamical systems involved in the relation
   *  \param time current time value
   */
  virtual void computeJacobianhOver_state(const siconos::algebra::BlockVector &state);

  /** Set a constant \f$ \nabla_{\lambda} g(\lambda) \f$ matrix for the system
   *
   *  \param newValue \f$ \nabla_{\lambda} g(\lambda) \f$ matrix
   *
   */
  void setConstantJacobiangOver_lambda(Eigen::Ref<siconos::algebra::SiconosMatrix> newValue);

  /** set a user-defined function to compute \f$ \nabla_{\lambda} g(\lambda) \f$ \f$
   *
   *  \param fct the user-defined function (std::function, lambda ...)
   */
  void setComputeJacobiangOver_lambdaFunction(
      const siconos::modeling::func_prototypes::FunctionV_M &fct);

  /** Computes \f$ \nabla_{\lambda} g(x, t, \lambda) \f$
   *  \param q coordinates of the dynamical systems involved in the relation
   *  \param time current time value
   */
  virtual void computeJacobiangOver_lambda(
      const Eigen::Ref<const siconos::algebra::SiconosVector> &lambda);

  /** default function to compute y, using the data from the Interaction and DS
   *
   *  \param time current time (not used)
   *  \param inter Interaction using this Relation
   *  \param level not used
   */
  void computeOutput(double time, Interaction &inter, unsigned int level = 0) override;

  /** default function to compute r, using the data from the Interaction and DS
   *
   *  \param time current time (not used)
   *  \param inter Interaction using this Relation
   *  \param level not used
   */
  void computeInput(double time, Interaction &inter, unsigned int level = 0) override;

  void computeJach(double time, Interaction &inter) override;

  void computeJacg(double time, Interaction &inter) override;

  /**
     return true if the relation requires the computation of residu

     \return true if residu are required, false otherwise
   */
  bool requireResidu() override { return true; }
};
}  // namespace siconos::modeling

#endif
