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

/*! \file FirstOrderNonLinearDS.hpp
  \brief First Order Non Linear Dynamical Systems
*/

#ifndef FIRSTORDERNONLINEARDS_H
#define FIRSTORDERNONLINEARDS_H

#include "DynamicalSystem.hpp"
#include "FunctionTypes.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosMemory.hpp"
#include "SiconosVector.hpp"

namespace siconos::algebra {

// class SiconosMatrix;
// class SiconosMatrix;
// class SiconosVector;
// class SiconosMemory;
}  // namespace siconos::algebra

namespace siconos::modeling {
/**
    General First Order Non Linear Dynamical Systems

    This class defines and computes a generic n-dimensional
    dynamical system of the form :

    \f[

    M \dot x = f(x,t) + r, \quad x(t_0) = x_0
    \f]
    where

    - \f$ x \in R^{n} \f$ is the state.
    - \f$ M \in R^{n\times n} \f$
    - \f$ r \in R^{n} \f$  is the input due to the Non Smooth Interaction.
    - \f$ f : R^{n} \times R  \mapsto  R^{n} \f$  a vector field.

    By default, the DynamicalSystem is considered to be an Initial Value Problem (IVP)
    and the initial conditions are given by

    \f[
    x(t_0)=x_0
    \f]

    To define a Boundary Value Problem, a pointer on a BoundaryCondition must be set.

    The right-hand side and its jacobian are defined as

    \f[
    rhs &=& \dot x =  M^{-1}(f(x,t)+ r) \\
    jacobianRhsOver_x &=& \nabla_x rhs(x,t) = M^{-1}\nabla_x f(x,t)
    \f]

    The following operators can set by user-defined functions (lambda ...)

    -  \f$ f(x,t) \f$
    -  \f$ \nabla_x f(x,t) \f$
    -  \f$ M(t) \f$

 */
class FirstOrderNonLinearDS : public DynamicalSystem {
 private:
  /** plugin signature */
  typedef void (*FNLDSPtrfct)(double, unsigned int, const double*, double*, unsigned int,
                              double*);

 protected:
  ACCEPT_SERIALIZATION(FirstOrderNonLinearDS);

  /** M matrix of the system (as a view onto memory)*/
  std::shared_ptr<siconos::algebra::MapType> MMatrix_view_{nullptr};

  /** internal storage (optional) for the MMatrix matrix */
  std::unique_ptr<std::vector<double>> MMatrix_internal_storage_{nullptr};

  /** function wrapper used to compute MMatrix */
  siconos::modeling::func_prototypes::FunctionS_M computeMMatrix_{nullptr};

  /** true if MMatrix is required/set and constant */
  bool hasConstantMMatrix_{false};

  /** True if MMatrix is required */
  bool hasMMatrix_{false};

  /** external forces applied to the system */
  std::shared_ptr<siconos::algebra::MapVectorType> fVector_view_{nullptr};

  /** internal (optional) storage used for external forces */
  std::unique_ptr<std::vector<double>> fVector_internal_storage_{nullptr};

  /** function wrapper used to compute external forces */
  siconos::modeling::func_prototypes::FunctionVS_V computefVector_{nullptr};

  /** True if external forces are taken into account and constant */
  bool hasConstantfVector_{false};

  /** True if external forces are taken into account */
  bool hasfVector_{false};

  /** \f$ \nabla_x f(x,t) \f$ (as a view onto memory)*/
  std::shared_ptr<siconos::algebra::MapType> jacobianfOver_x_view_{nullptr};

  /** internal storage (optional) for \f$ \nabla_x f(x,t) \f$ matrix */
  std::unique_ptr<std::vector<double>> jacobianfOver_x_internal_storage_{nullptr};

  /** function wrapper used to compute \f$ \nabla_x f(x,t) \f$ */
  siconos::modeling::func_prototypes::FunctionVS_M computejacobianfOver_x_{nullptr};

  /** true if \f$ \nabla_x f(x,t) \f$ is required/set and constant */
  bool hasConstantJacobianfOver_x_{false};

  /** True if \f$ \nabla_x f(x,t) \f$ is required */
  bool hasJacobianfOver_x_{false};

  /** A buffer to store f at previous time step \f$ (f(x_k,t_k) \f$) */
  std::shared_ptr<siconos::algebra::SiconosVector> fbuffer_{nullptr};

  /** Gradient of \f$ f(x,t) \f$ with respect to \f$ x \f$ */
  std::shared_ptr<siconos::algebra::SiconosMatrix> jacobianfOver_x_{nullptr};

  /**  the previous r vectors */
  siconos::algebra::SiconosMemory rMemory_;

  /**
      Copy of M Matrix, LU-factorized, used to solve systems like Mx = b with LU-factorization.
      (Warning: may not exist, used if we need to avoid factorization in place of M) */
  std::shared_ptr<siconos::algebra::SiconosDenseLUMatrix> LU_M_{nullptr};

  bool hasLU_M_{false};

  /** True if ALL components (f, \f$ \nablaf\f$, ...) are time independant or not used. Any
   * call to setCompute... turns this to false */
  bool isTimeInvariant_{true};

  /** flag to handle rhs init at first call for linear systems*/
  bool isFirstCall_{true};

  ///** default constructor */
  // FirstOrderNonLinearDS() = default;

 public:
  /** constructor from initial state, leads to \f$ \dot x = r \f$
   *
   *  \param newX0 initial state
   */
  FirstOrderNonLinearDS(Eigen::Ref<siconos::algebra::SiconosVector> x0);

  /** Copy constructor
   *
   *  \param FONLDS the FirstOrderNonLinearDS to copy
   */
  FirstOrderNonLinearDS(const FirstOrderNonLinearDS& FONLDS);

  /** destructor */
  virtual ~FirstOrderNonLinearDS() noexcept = default;

  /** allocate (if needed)  and compute rhs and its jacobian.
   *
   *  \param time of initialization
   */
  virtual void initRhs(double time) override;

  /** set nonsmooth input to zero
   *
   *  \param level input-level to be initialized.
   */
  void initializeNonSmoothInput(unsigned int level) override;

  /*  \return a read-only view on the matrix M(t) */
  inline auto MMatrix() const {
    return siconos::algebra::ConstMapType(MMatrix_view_->data(), MMatrix_view_->rows(),
                                          MMatrix_view_->cols());
  }

  /** Set a constant matrix M for the system
   *
   *  \param newValue MMatrix matrix
   */
  void setConstantMMatrix(Eigen::Ref<siconos::algebra::SiconosDenseMatrix> newValue);

  /** \return True if matrix M(t) has been set (i.e. different from identity) */
  bool hasMMatrix() const { return hasMMatrix_; }

  /** set a user-defined function to compute the matrix M(t)
   *
   *  \param fct the user-defined function (std::function, lambda ...)
   */
  void setComputeMMatrixFunction(const siconos::modeling::func_prototypes::FunctionS_M& fct);

  /** to compute the matrix M(t)
   *
   *  \param time current time value
   */
  virtual void computeMMatrix(double time);

  /** \return  a read-only view on f(x,t)  */
  inline auto fVector() const {
    return siconos::algebra::ConstMapVectorType(fVector_view_->data(), fVector_view_->size());
  }

  /** \return  a read-only view on fbuffer (value of f(x,t) at previous time step)  */
  inline auto fold() const {
    return siconos::algebra::ConstMapVectorType(fbuffer_->data(), fbuffer_->size());
  }

  /** set a constant f(x,t)  vector
   *
   *  \param newv f(x,t) vector
   */
  void setConstantfVector(Eigen::Ref<siconos::algebra::SiconosVector> newv);

  /** True if f(x,t) is taken into account */
  bool hasfVector() const { return hasfVector_; }

  /** set a user-defined function to compute f(x,t)
   *  \param fct the user-defined function (std::function, lambda ...)
   */
  void setComputefVectorFunction(const siconos::modeling::func_prototypes::FunctionVS_V& fct);

  /** Update  f(x,t)
   *
   * \param state state vector
   *  \param time the current time
   */
  void computefVector(const Eigen::Ref<siconos::algebra::SiconosVector>& state, double time);

  /** \return a read-only view on \f$ \nabla_xf(x,t) \f$ matrix */
  inline auto jacobianfOver_x() const {
    return siconos::algebra::ConstMapType(jacobianfOver_x_view_->data(),
                                          jacobianfOver_x_view_->rows(),
                                          jacobianfOver_x_view_->cols());
  }

  /** Set a constant \f$ \nabla_xf(x,t) \f$ matrix
   *
   *  \param newValue jacobianfOver_x matrix
   *
   */
  void setConstantJacobianfOver_x(Eigen::Ref<siconos::algebra::SiconosDenseMatrix> newValue);

  /** \return True if \f$ \nabla_xf(x,t) \f$ matrix has been set */
  bool hasJacobianfOver_x() const { return hasJacobianfOver_x_; }

  /** set a user-defined function to compute \f$ \nabla_xf(x,t) \f$
   *
   *  \param fct the user-defined function (std::function, lambda ...)
   */
  void setComputeJacobianfOver_xFunction(
      const siconos::modeling::func_prototypes::FunctionVS_M& fct);

  /** to compute  \f$ \nabla_xf(x,t) \f$
   *  \param state x vector
   *  \param time the current time
   */
  void computeJacobianfOver_x(const Eigen::Ref<siconos::algebra::SiconosVector>& state,
                              double time);

  /** reset the state to the initial state */
  void resetToInitialState() override;

  /** Change initial state (x0)
   * \param newX0 the new initial state
   *  warning: shared memory (eigen view) between newX0 and x0 attribute
   */
  void resetInitialState(Eigen::Ref<siconos::algebra::SiconosVector> newX0);

  /** Reset initial state (x0) components to a constant given value
   * \param val input value such that x0(i) = val
   * warning: reset also state to x0.
   */
  void resetInitialStateToValue(double val);

  /** \return a shared pointer to the initial state.
   *  This should be used with care (or not used at all ...)
   *  since the returned object will share memory with the internal value of x0 ...
   */
  auto x0_ptr() { return x0_view_; }
  // FP:  Used only in MatrixIntegrator for ZOH and friends. This should be avoided.

  /** update right-hand side for the current state
   *
   *  \param time of interest
   */
  virtual void computeRhs(double time) override;

  /** update  \f$ \nabla_x rhs \f$  for the current state
   *
   *  \param time of interest
   */
  void computeJacobianRhsOver_x(double time) override;

  /** reset non-smooth part of the rhs (i.e. r), for all 'levels' */
  virtual void resetAllNonSmoothParts() override;

  /** set nonsmooth part of the rhs (i.e. r) to zero for a given level
   *
   *  \param level
   */
  virtual void resetNonSmoothPart(unsigned int level) override;

  /** \return LU-factorization of the M matrix (pointer link) */
  inline auto LU_M() const { return LU_M_; }

  /** \return True if LU factorization of M is available */
  bool hasLU_M() const { return hasLU_M_; }

  /** get all the values of the state vector r stored in memory
   *
   *  \return a memory vector
   */
  inline const siconos::algebra::SiconosMemory& rMemory() const { return rMemory_; }

  /**
     initialize the SiconosMemory objects: reserve memory for i
      vectors in memory and reset all to zero.

      \param steps the size of the SiconosMemory (i)
  */
  void initMemory(unsigned int steps) override;

  /** push the current values of x and r in memory (index 0 of memory is the last inserted
   *  vector) xMemory and rMemory,
   */
  void swapInMemory() override;

  /**
     Call all plugged-function to initialize plugged-object values

     \param time value
  */
  virtual void updatePlugins(double time) override;

  /** print the data of the dynamical system on the standard output
   */
  void display(bool brief = true) const override;

  /** \return True if ALL components (f, \f$ \nablaf\f$, ...) are time independant or not used.
   * Any call to setCompute... turns this to false */
  auto isTimeInvariant() const { return isTimeInvariant_; }

  Type acceptType(types::FindType& ft) const override { return ft.visit(*this); }
};
}  // namespace siconos::modeling
#endif
