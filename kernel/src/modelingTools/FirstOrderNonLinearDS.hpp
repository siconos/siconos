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

    M \dot x = f(x,t,z) + r, \quad x(t_0) = x_0
    \f]
    where

    - \f$ x \in R^{n} \f$ is the state.
    - \f$ M \in R^{n\times n} \f$  a "mass matrix"
    - \f$ r \in R^{n} \f$  the input due to the Non Smooth Interaction.
    - \f$ z \in R^{zSize} \f$  is a vector of arbitrary algebraic
    variables, some sort of discret state.  For example, z may be used
    to set some perturbation parameters, to control the system (z
    set by actuators) and so on.

    - \f$ f : R^{n} \times R  \mapsto  R^{n} \f$  the vector field.

    By default, the DynamicalSystem is considered to be an Initial Value Problem (IVP)
    and the initial conditions are given by

    \f[
    x(t_0)=x_0
    \f]

    To define a Boundary Value Problem, a pointer on a BoundaryCondition must be set.

    The right-hand side and its jacobian (from base class) are defined as

    \f[
    rhs &=& \dot x =  M^{-1}(f(x,t,z)+ r) \\
    jacobianRhsx &=& \nabla_x rhs(x,t,z) = M^{-1}\nabla_x f(x,t,z)
    \f]

    The following operators can be plugged, in the usual way (see User Guide)

    -  \f$ f(x,t,z) \f$
    -  \f$ \nabla_x f(x,t,z) \f$
    -  \f$ M(t) \f$

 */
class FirstOrderNonLinearDS : public DynamicalSystem {
 public:
  using FunctionVS_M =
      std::function<void(const Eigen::Ref<const siconos::algebra::SiconosVector> &,
                         double,
                         Eigen::Ref<siconos::algebra::MapType>)>;

 private:
  /** plugin signature */
  typedef void (*FNLDSPtrfct)(double, unsigned int, const double *, double *, unsigned int,
                              double *);

 protected:
  ACCEPT_SERIALIZATION(FirstOrderNonLinearDS);

  /** Matrix coefficient of \f$ \dot x \f$ */
  std::shared_ptr<siconos::algebra::SiconosMatrix> _M{nullptr};

  /** value of f(x,t,z) */
  std::shared_ptr<siconos::algebra::SiconosVector> _f{nullptr};

  /** to store f(x_k,t_k,z_k)*/
  std::shared_ptr<siconos::algebra::SiconosVector> _fold{nullptr};

  /** view of \f$  \nabla_x f: (x,t) \f$ */
  std::shared_ptr<siconos::algebra::MapType> jacobianfx_view_{nullptr};

  /** internal storage for the jacobian view*/
  std::unique_ptr<std::vector<double>> jacobianfx_internal_storage_{nullptr};

  /** function wrapper used to compute jacobianfx */
  FunctionVS_M computejacobianfx_{nullptr};

  /** true if jacobianfx is constant */
  bool hasConstantJacobianfx_{false};

  /** True if jacobianfx is required */
  bool hasJacobianfx_{false};

  /** DynamicalSystem plug-in to compute f(x,t,z)
   *
   *  \param current time
   *  \param size of the vector _x
   *  \param[in,out] pointer to the first element of the vector _x
   *  \param[in,out] the pointer to the first element of the vector _f
   *  \param the size of the vector _z
   *  \param a vector of parameters _z
   */
  std::shared_ptr<siconos::plugins::PluggedObject> _pluginf{nullptr};

  std::shared_ptr<siconos::plugins::PluggedObject> _pluginM{nullptr};

  /**  the previous r vectors */
  siconos::algebra::SiconosMemory _rMemory;

  /**
      Copy of M Matrix, LU-factorized, used to solve systems like Mx = b with LU-factorization.
      (Warning: may not exist, used if we need to avoid factorization in place of M) */
  std::shared_ptr<siconos::algebra::SiconosMatrix> _invM{nullptr};

  /** default constructor */
  FirstOrderNonLinearDS() = default;

  /** Reset the PluggedObjects */
  void _zeroPlugin() override;

 public:
  // ===== CONSTRUCTORS =====

  /** constructor from initial state, leads to \f$ \dot x = r \f$
   *
   *  \param newX0 initial state
   *
   *  \warning you need to set explicitely the plugin for f and its jacobian if needed (e.g. if
   *  used with an EventDriven scheme)
   */
  FirstOrderNonLinearDS(std::shared_ptr<siconos::algebra::SiconosVector> newX0);

  /** constructor from initial state, leads to \f$ \dot x = r \f$
   *
   *  \param newX0 initial state
   *
   *  \warning you need to set explicitely the plugin for f and its jacobian if needed (e.g. if
   *  used with an EventDriven scheme)
   */
  FirstOrderNonLinearDS(Eigen::Ref<siconos::algebra::SiconosVector> &newX0);

  // /** constructor from initial state and f (plugins),  \f$ \dot x = f(x, t, z) + r \f$
  //  *
  //  *  \param newX0 initial state
  //  *  \param fPlugin name of the plugin function to be used for f(x,t,z)
  //  *  \param jacobianfxPlugin name of the plugin to be used for the jacobian of f(x,t,z)
  //  */
  // FirstOrderNonLinearDS(std::shared_ptr<siconos::algebra::SiconosVector> newX0,
  //                       const std::string &fPlugin, const std::string &jacobianfxPlugin);

  /** Copy constructor
   *
   *  \param FONLDS the FirstOrderNonLinearDS to copy
   */
  FirstOrderNonLinearDS(const FirstOrderNonLinearDS &FONLDS);

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

  /** reset the state to the initial state */
  void resetToInitialState() override;

  /** update right-hand side for the current state
   *
   *  \param time of interest
   */
  virtual void computeRhs(double time) override;

  /** update  \f$ \nabla_x rhs \f$  for the current state
   *
   *  \param time of interest
   */
  void computeJacobianRhsx(double time) override;

  /** reset non-smooth part of the rhs (i.e. r), for all 'levels' */
  virtual void resetAllNonSmoothParts() override;

  /** set nonsmooth part of the rhs (i.e. r) to zero for a given level
   *
   *  \param level
   */
  virtual void resetNonSmoothPart(unsigned int level) override;

  /** returns a pointer to M, matrix coeff. on left-hand side
   */
  inline std::shared_ptr<siconos::algebra::SiconosMatrix> M() const { return _M; }

  /** set M, matrix coeff of left-hand side (pointer link)
   *
   *  \param newM the new M matrix
   */
  inline void setMPtr(std::shared_ptr<siconos::algebra::SiconosMatrix> newM) { _M = newM; }

  /** get the inverse of LU fact. of M operator (pointer link)
   *
   *  \return pointer to a SiconosMatrix
   */
  inline std::shared_ptr<siconos::algebra::SiconosMatrix> invM() const { return _invM; }

  /** returns f(x,t,z) (pointer link)
   */
  inline std::shared_ptr<siconos::algebra::SiconosVector> f() const { return _f; }

  /** set f(x,t,z) (pointer link)
   *
   *  \param newPtr a std::shared_ptr<siconos::algebra::SiconosVector>
   */
  inline void setFPtr(std::shared_ptr<siconos::algebra::SiconosVector> newPtr) { _f = newPtr; }

  /*  \return \f$  \nabla_x f: (x,t) \in R^{n} \times R \mapsto R^{n \times n} \f$ */
  inline virtual std::shared_ptr<siconos::algebra::MapType> jacobianfx() const {
    return jacobianfx_view_;
  }

  /** \return a view onto \f$  \nabla_x f: (x,t) \f$  */
  inline siconos::algebra::MapType &jacobianfx_view() const { return *jacobianfx_view_; }

  /** set  \f$  \nabla_x f: (x,t) \f$ as a constant matrix.
   *
   *  \param newValue new jacobian matrix
   *
   */
  void setConstantJacobianfx(Eigen::Ref<siconos::algebra::SiconosMatrix> newValue);

  /** \return True if jacobian matrix is taken into account in the system */
  bool hasJacobianfx() const { return hasJacobianfx_; }

  /** set a user-defined function to compute \f$  \nabla_x f: (x,t) \f$
   *
   *  \param fct the user-defined function (std::function, lambda ...)
   */
  void setComputeJacobianfxFunction(const FunctionVS_M& fext_func);

  /** default function to compute \f$  \nabla_x f: (x,t) \f$
   *
   *  \param time the current time
   *  \param position the current state vector
   */
  void computeJacobianfx(double time, Eigen::Ref<siconos::algebra::SiconosVector> position);

  /** Default function to compute  \f$  \nabla_x f: (x,t) \in R^{n}
   *  \times R \mapsto R^{n \times n} \f$ with x different from
   *  current saved state.
   *
   *  \param time instant used in the computations
   *  \param state a siconos::algebra::SiconosVector to store the resuting value
   */

  /** get all the values of the state vector r stored in memory
   *
   *  \return a memory vector
   */
  inline const siconos::algebra::SiconosMemory &rMemory() const { return _rMemory; }

  /** returns previous value of rhs -->OSI Related!!*/
  inline std::shared_ptr<siconos::algebra::SiconosVector> fold() const { return _fold; }

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

  // --- setters for functions to compute plugins ---

  /** to set a specified function to compute M
   *
   *  \param pluginPath the complete path to the plugin
   *  \param functionName function name to use in this library
   */
  void setComputeMFunction(const std::string &pluginPath, const std::string &functionName);

  /** set a specified function to compute M
   *
   *  \param fct a pointer on the plugin function
   */
  void setComputeMFunction(siconos::plugins::FPtr1 fct);

  /** to set a specified function to compute f(x,t)
   *
   *  \param pluginPath the complete path to the plugin
   *  \param functionName the function name to use in this library
   */
  void setComputeFFunction(const std::string &pluginPath, const std::string &functionName);

  /** set a specified function to compute the vector f
   *
   *  \param fct a pointer on the plugin function
   */
  void setComputeFFunction(siconos::plugins::FPtr1 fct);

  // --- compute plugin functions ---

  /** Default function to compute  \f$  M: (x,t) \f$
   *
   *  \param time time instant used in the computations
   */
  void computeM(double time);

  /** Default function to compute \f$ f: (x,t) \f$
   *
   *  \param time time instant used in the computations
   */
  // virtual void computef(double time);

  /** function to compute \f$ f: (x,t) \f$
   *
   *  \param time time instant used in the computations
   *  \param state x value
   */
  virtual void computef(double time, std::shared_ptr<siconos::algebra::SiconosVector> state);

  /** Get _pluginf
   *
   *  \return a std::shared_ptr<siconos::plugins::PluggedObject>
   */
  inline std::shared_ptr<siconos::plugins::PluggedObject> getPluginF() const {
    return _pluginf;
  };

  /** Get _pluginM
   *
   *  \return a std::shared_ptr<siconos::plugins::PluggedObject>
   */
  inline std::shared_ptr<siconos::plugins::PluggedObject> getPluginM() const {
    return _pluginM;
  };

  /** print the data of the dynamical system on the standard output
   */
  void display(bool brief = true) const override;

  // visitors hook
  void acceptSP(std::shared_ptr<siconos::internal::SiconosVisitor> tourist) const override;
  Type acceptType(types::FindType &ft) const override { return ft.visit(*this); }
};
}  // namespace siconos::modeling
#endif
