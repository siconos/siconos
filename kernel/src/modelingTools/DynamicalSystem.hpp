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
 * Unless required by applicable law or agreed to in writing, softwa∏re
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
/*! \file DynamicalSystem.hpp
  Abstract class - General interface for all Dynamical Systems.
*/

#ifndef DynamicalSystem_H
#define DynamicalSystem_H

#include <memory>
#include <vector>

#include "DynamicalSystemVisitor.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosMemory.hpp"
#include "SiconosVector.hpp"
#include "StorageTools.hpp"
#include "TypeName.hpp"  // visitor to get ds type

namespace siconos::modeling {

/**

    Abstract interface to Dynamical Systems

    This class is used to describe dynamical systems of the form :

    \f$ g(\dot x, x, t) = 0 \f$

    where

    - \f$ x \in R^{n} \f$ is the state.
    - \f$ g : R^{n} \times R  \to  R^{n}   \f$ .

    By default, the DynamicalSystem is considered to be an Initial Value
    Problem (IVP) and the initial conditions are given by

    \f$ x(t_0)=x_0 \f$

    Under some specific conditions, the system can be written as:

    \f$ \dot x = rhs(x, t) \f$

    In that case, \f$ \nabla_{\dot x} g \f$ must be invertible.

*/
class DynamicalSystem {
 public:
  /** List of indices used to save tmp work vectors
   * The last value is the size of the present list, so you HAVE to leave it at
   * the end position.
   */
  enum class DSWorkVectorId {
    local_buffer,
    freeresidu,
    free,
    acce_memory,
    acce_like,
    sizeWorkV
  };

 private:
  ACCEPT_SERIALIZATION(DynamicalSystem);

  /** used to set ds number */
  static size_t __count;

 protected:
  /** An id number for the DynamicalSystem */
  size_t number_{__count++};

  /** the dimension of the system (\e ie size of the state vector x) */
  siconos::algebra::Index x_size_{0};

  /** initial state of the system */
  siconos::algebra::DenseVectorStorage x0_storage_{std::monostate{}};

  template <typename F>
  decltype(auto) use_x0(F&& f) const {
    return siconos::algebra::visitStorage(x0_storage_, std::forward<F>(f), "x0_storage_");
  }
  template <typename F>
  decltype(auto) use_x0(F&& f) {
    return siconos::algebra::visitStorage(x0_storage_, std::forward<F>(f), "x0_storage_");
  }

  /** the input vector due to the non-smooth law \f$ r \in R^{n} \f$
   * (multiplier, force, ...)
   * \remark V.A. 17/09/2011 :
   * This should be a VectorOfVectors as for state_x_ when higher relative degree
   * systems will be simulated
   */
  std::shared_ptr<siconos::algebra::SiconosVector> rVector_{nullptr};

  /** state of the system,
   *  \f$  x \in R^{n} \f$ - With state_x_[0]= \f$ x \f$ , state_x_[1]= \f$ \dot{x} \f$ . */
  std::vector<std::shared_ptr<siconos::algebra::SiconosVector>> state_x_ = {nullptr, nullptr};

  /** jacobian according to x of the right-hand side (\f$ rhs = \dot x = M^{-1}
      f(x,t) + r \f$) */
  siconos::algebra::SiconosVector jacobianRhsOver_x_;
  // Note FP: used only for EventDriven with Lsodar or similar integrators, based
  // on fortran wrappers. We only need to get this jacobian as a vector
  // which contains the columns of the real jacobian.

  /** Internal buffer, used in rhs or jacobian rhs computation. Do not trust its content!
  Empty by default
   */
  siconos::algebra::SiconosVector buffer_;

  /** False if \f$ rhs = \dot x = M^{-1} f(x,t) + r \f$) needs to be updated */
  bool is_jacobianRhsOver_x_uptodate_{false};

  /** the  previous state vectors stored in memory */
  siconos::algebra::SiconosMemory xMemory_;

  /** number of previous states stored in memory */
  unsigned int stepsInMemory_{1};

  // ===== CONSTRUCTORS =====

  /** default constructor */
  DynamicalSystem() = default;

  /** minimal constructor, from state dimension
      result in \f$ \dot x = r \f$
   *  \param dimension size of the system (n)
   */
  DynamicalSystem(siconos::algebra::Index dimension) : x_size_(dimension) {};

  /** Copy constructor
   * \param ds the DynamicalSystem to copy
   */
  DynamicalSystem(const DynamicalSystem& ds);

  // Rule of five
  DynamicalSystem& operator=(const DynamicalSystem&) = delete;
  DynamicalSystem(DynamicalSystem&&) = delete;
  DynamicalSystem& operator=(DynamicalSystem&&) = delete;

 public:
  /** destructor */
  virtual ~DynamicalSystem() noexcept = default;

  /** \return a read-only view on the initial state vector */
  /**  \return a read-only reference on the initial state */
  Eigen::Ref<const siconos::algebra::SiconosVector> x0() const {
    return use_x0(
        [](auto const& v) { return Eigen::Ref<const siconos::algebra::SiconosVector>(v); });
  }
  /** True if x0 is defined */
  bool hasX0() const { return !std::holds_alternative<std::monostate>(x0_storage_); }

  /** allocate (if needed)  and compute rhs and its jacobian.
   *
   *  \param time of initialization
   */
  virtual void initRhs(double time) = 0;

  /** set nonsmooth input to zero
   *
   *  \param level input-level to be initialized.
   */
  virtual void initializeNonSmoothInput(unsigned int level) = 0;

  /** compute all component of the dynamical system, for the current state.
   *
   *  \param time current time (the one used to update ds component)
   */
  void update(double time);

  /** update right-hand side for the current state
   *
   *  \param time of interest
   */
  virtual void computeRhs(double time) = 0;

  /** update \f$ \nabla_x rhs \f$ for the current state
   *
   *  \param time of interest
   */
  virtual void computeJacobianRhsOver_x(double time) = 0;

  /** reset nonsmooth part of the rhs, for all 'levels' */
  virtual void resetAllNonSmoothParts() = 0;

  /** set nonsmooth part of the rhs to zero for a given level
   *
   *  \param level
   */
  virtual void resetNonSmoothPart(unsigned int level) = 0;

  /** returns the id of the dynamical system */
  inline size_t number() const { return number_; }

  /** set the id of the DynamicalSystem
   *
   *  \return the previous value of number
   */
  inline size_t setNumber(size_t new_number) {
    size_t old_n = number_;
    number_ = new_number;
    return old_n;
  }

  /** returns the size of the vector state x */
  inline auto x_size() const { return x_size_; }

  /** Fake function to access the dimension used to allocate iteration matrix
   *  in the integrators
   * Usually equal to x_size or ndof but might be larger (see fem).
   * \return number of degrees of freedom of the system
   */
  inline virtual siconos::algebra::Index real_size() const { return x_size_; }

  /**
      returns the dimension of the system
      (depends on system type, e.g. n for first order,
      ndof for Lagrangian).
   */
  virtual inline siconos::algebra::Index dimension() const { return x_size_; };

  /** \return a read-only view on the state vector (size=dimension()) */
  inline const siconos::algebra::ConstMapVectorType x_read() const {
    return siconos::algebra::ConstMapVectorType(state_x_[0]->data(), state_x_[0]->size());
  }

  /** \return the state vector \f$ x \f$ */
  inline std::shared_ptr<siconos::algebra::SiconosVector> x() const { return state_x_[0]; }

  // /** generalized coordinates of the system (vector of size dimension())
  //  *
  //  *  \return pointer on a siconos::algebra::SiconosVector
  //  */
  inline siconos::algebra::SiconosVector& x_python() const { return *(state_x_[0]); }

  /** \return r vector (input due to nonsmooth behavior) */
  inline std::shared_ptr<siconos::algebra::SiconosVector> r() const { return rVector_; }

  /** \return r vector (input due to nonsmooth behavior) */
  inline siconos::algebra::SiconosVector& r_python() const { return *rVector_; }

  /** \return the right-hand side vector (i.e. \f$ \dot x \f$) */
  inline std::shared_ptr<siconos::algebra::SiconosVector> rhs() const { return state_x_[1]; }

  /** \return a read-only view on the right-hand side vector (i.e. \f$ \dot x \f$) */
  inline const siconos::algebra::ConstMapVectorType rhs_read() const {
    return siconos::algebra::ConstMapVectorType(state_x_[1]->data(), state_x_[1]->size());
  }

  /** \return a read-only access to \f$ \nabla_x rhs()\f$, as a vector =[col0, col1, ...]
   */
  const siconos::algebra::SiconosVector& jacobianRhsOver_x() const {
    return jacobianRhsOver_x_;
  }

  /** \return a writable vector which represents \f$ \nabla_x rhs()\f$, as a vector =[col0,
   * col1, ...]
   */
  siconos::algebra::SiconosVector& jacobianRhsOver_x_readwrite() { return jacobianRhsOver_x_; }

  /** get all the values of the state vector x stored in a SiconosMemory object
   *  (not const due to LinearSMC::actuate)
   *
   *  \return a reference to the SiconosMemory object
   */
  inline siconos::algebra::SiconosMemory& xMemory() { return xMemory_; }

  /** get all the values of the state vector x stored in a SiconosMemory object
   *
   *  \return a const reference to the SiconosMemory object
   */
  inline const siconos::algebra::SiconosMemory& xMemory() const { return xMemory_; }

  /** returns the number of step saved in memory for state vector
   *
   *  \return int
   */
  inline auto stepsInMemory() const { return stepsInMemory_; }

  /** set number of steps to be saved
   *
   *  \param steps
   */
  inline void setStepsInMemory(unsigned int steps) { stepsInMemory_ = steps; }

  /** initialize the SiconosMemory objects: reserve memory for i vectors in
   * memory and reset all to zero.
   *
   *  \param steps the size of the SiconosMemory (i)
   */
  virtual void initMemory(unsigned int steps);

  /** push the current values of x and r in memory (index 0 of memory is the
   * last inserted vector) xMemory and rMemory,
   */
  virtual void swapInMemory() = 0;

  /** call all plugged functions for the current state
   *
   *  \param time  the current time
   */
  virtual void updatePlugins(double time) {};

  /** reset the global DynamicSystem counter (for ids)
   *
   *  \return the previous value of count
   */
  static inline size_t resetCount(size_t new_count = 0) {
    size_t old_count = __count;
    __count = new_count;
    return old_count;
  };

  /** @brief Reset initial state.
   *
   * Warning : deep copy of the provided vector into internal attribute
   *
   * @param newValue external vector to be copied. Its size must match dimension()
   * @param tag Pass siconos::algebra::copy_t to select this overload (rather than alias
  version)
   *
   */
  void setX0(const siconos::algebra::SiconosVector& newValue, siconos::algebra::CopyTag tag);

  /** @brief set a constant external forces vector
   *
   * Warning : This method does NOT copy the data. Instead, it creates an Eigen::Map
   * pointing directly to the memory provided by the argument.
   *
   * This means:
   *  - ownership stays external
   *  - modifications to the original vector are reflected inside the class
   *
   * @param newValue external vector to be copied. Its size must match dimension()
   * @param tag Pass siconos::algebra::alias_t to select this overload
   *        (rather than copy version)
   *
   */
  void setX0(Eigen::Ref<siconos::algebra::SiconosVector> newValue,
             siconos::algebra::AliasTag tag);

  /** reset the state x() to the initial state x0 */
  virtual void resetToInitialState();

  /** \return true if the system is linear
   */
  virtual bool isLinear() const { return false; };

  /** print the data of the dynamical system on the standard output
   */
  virtual void display(bool brief = true) const = 0;

  // accept function used to call visit(shared_ptr<DS>)
  virtual void acceptSP(dynamical_systems::Visitor&) {
    throw std::logic_error(
        "this class derived from DynamicalSystem does not accept a visitor for shared "
        "pointers");
  };

  // accept function used to call visit(DS)
  virtual void accept(dynamical_systems::Visitor&) const {
    throw std::logic_error("accept: no visitor defined");
  };

  virtual Type acceptType(siconos::types::FindType& ft) const = 0;
};
}  // namespace siconos::modeling
#endif  // DynamicalSystem_H
