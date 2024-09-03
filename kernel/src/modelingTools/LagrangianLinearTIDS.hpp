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

/*! \file LagrangianLinearTIDS.hpp */
#ifndef LAGRANGIANTIDS_H
#define LAGRANGIANTIDS_H

#include "LagrangianDS.hpp"
#include "matrix_wrapper.h"
#include "vector_wrapper.h"
#include "SiconosPointers.hpp"

namespace siconos::modeling {

/**
    Lagrangian Linear Systems with time invariant coefficients -

    \f$ M\dot v + Cv + Kq = F_{ext}(t,z) + p \f$

    The class LagrangianLinearTIDS  allows to define  and compute a generic ndof-dimensional
   Lagrangian Linear Time Invariant Dynamical System of the form:

    \f[
    M \ddot q + C \dot q + K q =  F_{ext}(t,z) + p,
    \f]

    where
    -  \f$ q \in R^{ndof} \f$ is the set of the generalized coordinates,
    - \f$ \dot q  \in R^{ndof} \f$  the velocity, i. e. the time derivative of
    the  generalized coordinates.
    - \f$ \ddot q  \in R^{ndof} \f$  the acceleration, i. e. the second time
    derivative of the  generalized coordinates.
    - \f$ p \in R^{ndof} \f$  the forces due to the Non Smooth Interaction. In
    particular case of Non Smooth evolution, the variable p contains the impulse
    and not the force.
    -  \f$ M \in  R^{ndof \times ndof} \f$ is the Mass matrix (access : mass()
    method).
    -  \f$ K \in  R^{ndof \times ndof} \f$ is the stiffness matrix (access : K()
    method).
    -  \f$ C \in  R^{ndof \times ndof} \f$ is the viscosity matrix (access : C()
    method).
    -  \f$ z \in R^{zSize} \f$  is a vector of arbitrary algebraic variables,
   some sort of discret state.

    The equation of motion is also shortly denoted as:

    \f[

    M(q,z) \dot v = F(v, q, t, z) + p

    \f]

    where
    -  \f$ F(v, q, t, z) \in R^{ndof} \f$ collects the total forces
    acting on the system, that is
    \f$ F(v, q, t, z) =  F_{ext}(t, z) -  Cv - Kq \f$.

    This vector is saved and may be accessed using forces() method.

    If required (e.g. for Event-Driven like simulation), reformulation as a
    first-order system is also available, and writes:

    - \f$ n= 2 ndof \f$
    - \f$ x = \left[\begin{array}{c}q \\ \dot q\end{array}\right] \f$
    - rhs given by:

    \f[

        rhs(x,t,z) = \left[\begin{array}{c}
        \dot q  \\
        \ddot q = M^{-1}\left[F_{ext}(t, z) - C \dot q - K q  + p \right] \\
        \end{array}\right]

    \f]

    Its jacobian is:

    \f[

      \nabla_{x}rhs(x,t) = \left[\begin{array}{cc}
      0   & I \\
      -M^{-1}K & -M^{-1}C \\
      \end{array}\right]

    \f]

    with the input due to the non smooth law:

    \f[
      r = \left[\begin{array}{c}0 \\ p \end{array}\right]

    \f]
*/
class LagrangianLinearTIDS : public LagrangianDS {
 protected:
  ACCEPT_SERIALIZATION(LagrangianLinearTIDS);

  /** stiffness matrix */
  std::shared_ptr<Matrix> _K{nullptr};

  /** damping matrix */
  std::shared_ptr<Matrix> _C{nullptr};

  /** default constructor */
  LagrangianLinearTIDS() = default;

 public:
  std::shared_ptr<siconos::algebra::MapType> _massMap {nullptr};
  /** constructor from initial state and all matrix operators.
   *
   *  \param q0 initial coordinates
   *  \param v0 initial velocity
   *  \param M mass matrix
   *  \param K stiffness matrix
   *  \param C damping matrix
   */
  LagrangianLinearTIDS(std::shared_ptr<siconos::algebra::SiconosVector> q0,
                       std::shared_ptr<siconos::algebra::SiconosVector> v0,
                       std::shared_ptr<Matrix> M, std::shared_ptr<Matrix> K,
                       std::shared_ptr<Matrix> C);

  /** constructor from initial state and mass matrix only. Leads to \f$ M\dot v
   *  = F_{ext}(t,z) + p \f$ .
   *
   *  \param q0 initial coordinates
   *  \param v0 initial velocity
   *  \param M mass matrix
   */
  LagrangianLinearTIDS(std::shared_ptr<siconos::algebra::SiconosVector> q0,
                       std::shared_ptr<siconos::algebra::SiconosVector> v0,
                       std::shared_ptr<Matrix> M)
      : LagrangianDS(q0, v0, M){};

  // LagrangianLinearTIDS(VectorWrapper q0,
  //                      VectorWrapper v0,
  //                      MatrixWrapper m)
  //     : LagrangianLinearTIDS{q0.get_shared_ptr(), v0.get_shared_ptr(), m.get_shared_ptr()}{};
      

/** constructor from initial state and mass matrix only. Leads to \f$ M\dot v
   *  = F_{ext}(t,z) + p \f$ .
   *
   *  \param q0 initial coordinates
   *  \param v0 initial velocity
   *  \param M mass matrix
   */
  LagrangianLinearTIDS(Eigen::Ref<siconos::algebra::SiconosVector> &q0,
                       Eigen::Ref<siconos::algebra::SiconosVector> &v0,
                       Eigen::Ref<siconos::algebra::SiconosMatrix> &M)
      : LagrangianDS(){
        // using MapType = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>>;
        this->_massMap = std::make_shared<siconos::algebra::MapType>(M.data(), M.rows(), M.cols());

        // Here there is copy, test purpose
        _mass = std::make_shared<siconos::algebra::SiconosMatrix>(M);
        // Solution qui marche (need to template DS with mass matrix size)
        // Eigen::Map<Eigen::Matrix<double, Rows, Cols>> tmp(NULL);
        // new (&(tmp)) Eigen::Map<Eigen::Matrix<double, Rows, Cols>>(M.data(), M.rows(), M.cols());
        // for (int i = 0; i < tmp.cols(); i++)
        // {
        //   tmp(0, i) = 403.;
        // }
        // tmp.display();
      };

  siconos::algebra::MapType& mass_python() {
    return *_massMap;
  }

  siconos::algebra::SiconosMatrix& mass2() {
    return *_mass;
  }

  /** destructor */
  ~LagrangianLinearTIDS() noexcept = default;

  /** allocate (if needed)  and compute rhs and its jacobian.
   *
   *  \param t time of initialization
   */
  void initRhs(double t) override;

  /** Compute \f$ F(v,q,t,z) \f$
   *
   *  \param time the current time
   *  \param q std::shared_ptr<siconos::algebra::SiconosVector>: pointers on q
   *  \param velocity std::shared_ptr<siconos::algebra::SiconosVector>: pointers
   * on velocity
   */
  void computeForces(double time, std::shared_ptr<siconos::algebra::SiconosVector> q,
                     std::shared_ptr<siconos::algebra::SiconosVector> velocity) override;

  /** \return the stiffness matrix (pointer link) */
  inline std::shared_ptr<Matrix> K() const { return _K; }

  /** set (copy) the value of the stiffness matrix
   *
   *  \param K new stiffness matrix
   */
  void setK(const Matrix &K);

  /** set stiffness matrix (pointer link)
   *
   *  \param newPtr pointer to the new Stiffness matrix
   */
  void setKPtr(std::shared_ptr<Matrix> newPtr);

  /** \return the damping matrix (pointer link) */
  inline std::shared_ptr<Matrix> C() const { return _C; }

  /** set (copy) the value of the damping matrix
   *
   *  \param C new damping matrix
   */
  void setC(const Matrix &C);

  /** set damping matrix (pointer link)
   *
   *  \param newPtr pointer to the new damping matrix
   */
  void setCPtr(std::shared_ptr<Matrix> newPtr);

  /** \return \f$ \nabla_qF(v,q,t,z) \f$ (pointer  link) */
  inline std::shared_ptr<Matrix> jacobianqForces() const override { return _K; }

  /** \return \f$ \nabla_{\dot q}F(v,q,t,z) \f$ (pointer  link) */
  inline std::shared_ptr<Matrix> jacobianvForces() const override { return _C; }

  /** \return true if the Dynamical system is linear.
   */
  virtual bool isLinear() override { return true; }

  /** print the data onto the screen
   */
  void display(bool brief = true) const override;

  Type acceptType(types::FindType &ft) const override { return ft.visit(*this); }
};
}  // namespace siconos::modeling
#endif  // LAGRANGIANTIDS_H
