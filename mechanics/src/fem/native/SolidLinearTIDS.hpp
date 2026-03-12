/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2026 INRIA.
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

/*! \file SolidLinearTIDS.hpp

 */
#ifndef SOLIDLAGRANGIANTIDS_H
#define SOLIDLAGRANGIANTIDS_H

#include "FiniteElementLinearTIDS.hpp"

/** Finite Element discretization of an elastic solids that inherits from FEM with time
 * invariant coefficients and exhibits sigma to enable easy plasticty implementation using a
 * stress relation.
 * - \f$M\dot v + K u +  B^T \sigma = F_{ext}(t,z) + p
 * v = \dot u
 * \dot \sigma = B v + \dot \epsilon_p \f$
 */
namespace siconos::mechanics::fem {
class SolidLinearTIDS : public FiniteElementLinearTIDS {
 protected:
  ACCEPT_SERIALIZATION(SolidLinearTIDS);

  /** Elasticity Matrix */
  siconos::algebra::SparseStorage elasticity_matrix_storage_{std::monostate{}};
  template <typename F>
  decltype(auto) use_elasticity_matrix(F&& f) const {
    return siconos::algebra::visitStorage(elasticity_matrix_storage_, std::forward<F>(f),
                                          "elasticity_matrix_storage_");
  }

  /** B matrix from FEM */
  siconos::algebra::SparseStorage B_matrix_storage_{std::monostate{}};
  template <typename F>
  decltype(auto) use_B_matrix(F&& f) const {
    return siconos::algebra::visitStorage(B_matrix_storage_, std::forward<F>(f),
                                          "B_matrix_storage_");
  }

  /** Initial stress of the system */
  std::unique_ptr<siconos::algebra::SiconosVector> sigma0_{nullptr};

  /** Initial strain of the system */
  std::unique_ptr<siconos::algebra::SiconosVector> epsilon_p0_{nullptr};

  /** stress of the system */
  std::unique_ptr<siconos::algebra::SiconosVector> sigma_{nullptr};

  siconos::algebra::blocks::SharedVector stress_{nullptr, nullptr, nullptr};

  siconos::algebra::blocks::SharedVector epsilon_p_ = {nullptr, nullptr, nullptr};
  std::shared_ptr<siconos::algebra::SiconosVector> plasticRate_ = {nullptr};

  /** total size of sigma tensor */
  siconos::algebra::Index dimStress_{0};

  /** memory of previous coordinates of the system */
  siconos::algebra::SiconosMemory stressMemory_;
  //  siconos::algebra::SiconosMemory stressRateMemory_;
  // siconos::algebra::SiconosMemory plasticDeformationMemory_;
  // siconos::algebra::SiconosMemory plasticRate_Memory_;

 public:
  /** @brief Standard constructor from FE mesh and material description
   *  @param[in] mesh the mesh that defined the spatial discretization
   *  @param[in] materials map describing the different material properties
   */
  SolidLinearTIDS(std::shared_ptr<Mesh> mesh, const std::map<int, const Material>& materials);

  /** destructor */
  ~SolidLinearTIDS() noexcept = default;

  /** @brief print system data to screen */
  void display(bool brief = true) const override;

  /** Fake function to access the dimension used to allocate iteration matrix
   *  in the integrators
   *  Here: ndof + size of stress
   * \return dimension + stress dimension
   */
  inline virtual siconos::algebra::Index real_size() const override {
    return dimension() + stressDimension();
  }

  /** @return  a read-only view on stress  */
  inline auto stress_read() const {
    return siconos::algebra::ConstMapVectorType(stress_[0]->data(), stress_[0]->size());
  }

  /** @return the stress (pointer link) */
  std::shared_ptr<siconos::algebra::SiconosVector> stress() const { return stress_[0]; }

  /** @return the stress (pointer link) */
  std::shared_ptr<siconos::algebra::SiconosVector> stressRate() const { return stress_[1]; }

  // To be used only for pybind11 !
  inline siconos::algebra::SiconosVector& stress_python() const { return *(stress_[0]); }

  /** @return  a read-only view on stress_rate  */
  inline auto stress_rate_read() const {
    return siconos::algebra::ConstMapVectorType(stress_[1]->data(), stress_[1]->size());
  }

  /** @return the stress_rate (pointer link) */
  std::shared_ptr<siconos::algebra::SiconosVector> stress_rate() const { return stress_[1]; }

  // To be used only for pybind11 !
  inline siconos::algebra::SiconosVector& stress_rate_python() const { return *(stress_[1]); }

  /** @return true if epsilon_p[level] is used (i.e. has been allocated)
   *  @param level expected level
   */
  bool has_epsilon_p(size_t level) const { return epsilon_p_[level] != nullptr; }

  /** @return  a read-only view on epsilon_p[level]
   * @param level required level
   */
  inline auto epsilon_p(size_t level) const {
    return siconos::algebra::ConstMapVectorType(epsilon_p_[level]->data(),
                                                epsilon_p_[level]->size());
  }

  /** @return  a read-only view on the plastic deformation */
  inline auto plasticDeformation_read() const {
    return siconos::algebra::ConstMapVectorType(epsilon_p_[0]->data(), epsilon_p_[0]->size());
  }

  /** @return the plastic deformation (pointer link) */
  std::shared_ptr<siconos::algebra::SiconosVector> plasticDeformation() const {
    return epsilon_p_[0];
  }

  /** @return  a read-only view on the plastic deformation rate */

  inline auto plasticRate_read() const {
    return siconos::algebra::ConstMapVectorType(epsilon_p_[1]->data(), epsilon_p_[1]->size());
  }

  /** @return the plastic deformation rate (pointer link) */
  std::shared_ptr<siconos::algebra::SiconosVector> plasticRate() const {
    return epsilon_p_[1];
  }

  /** @brief set the elasticity matrix.
   *
   * Warning : This method does NOT copy the data. Instead, it creates an Eigen::Map
   * pointing directly to the memory provided by the argument.
   *
   * This means:
   *  - ownership stays external
   *  - modifications to the original vector are reflected inside the class
   * @param S new elasticity matrix
   * @param tag Pass siconos::algebra::alias_t to select this overload
   *        (rather than copy version)
   */
  void setElasticityMatrix(Eigen::Map<siconos::algebra::SiconosSparseMatrix>& S,
                           siconos::algebra::AliasTag);

  /** @brief Set a constant elasticity matrix for the system
   *
   *  Warning : deep copy of the provided object into internal attribute
   *
   *  @param newValue elasticity matrix
   *  @param tag Pass siconos::algebra::copy_t to select this overload (rather than alias
   *
   */
  template <typename T>
  void setElasticityMatrix(T&& newValue, siconos::algebra::CopyTag) {
    static_assert(std::is_same_v<std::decay_t<T>, siconos::algebra::SiconosSparseMatrix>,
                  "Type must be SiconosSparseMatrix");

    if (newValue.rows() != ndof_ || newValue.cols() != ndof_)
      throw std::invalid_argument("Input elasticity matrix has wrong dimensions");

    elasticity_matrix_storage_ =
        std::make_unique<siconos::algebra::SiconosSparseMatrix>(std::forward<T>(newValue));
  }

  /** @return a read-only ref onto the elasticity matrix */
  inline Eigen::Ref<const siconos::algebra::SiconosSparseMatrix> elasticityMatrix() const {
    return use_elasticity_matrix(
        [](Eigen::Ref<const siconos::algebra::SiconosSparseMatrix> M) {
          return Eigen::Ref<const siconos::algebra::SiconosSparseMatrix>(M);
        });
  }

  /*  @return a reference to the elasticity matrix

    pybind11 use only!
  */
  const siconos::algebra::SiconosSparseMatrix& elasticityMatrix_py() const {
    return siconos::algebra::visitStorage<siconos::algebra::AccessMode::OwnedOnly>(
        elasticity_matrix_storage_,
        [](const auto& M) -> const siconos::algebra::SiconosSparseMatrix& { return M; },
        "elasticity_matrix_storage_");
  }

  /** True if elasticity matrix is defined */
  bool hasElasticityMatrix() const {
    return !std::holds_alternative<std::monostate>(elasticity_matrix_storage_);
  }

  /** @brief set the B matrix.
   *
   * Warning : This method does NOT copy the data. Instead, it creates an Eigen::Map
   * pointing directly to the memory provided by the argument.
   *
   * This means:
   *  - ownership stays external
   *  - modifications to the original vector are reflected inside the class
   * @param S new B matrix
   * @param tag Pass siconos::algebra::alias_t to select this overload
   *        (rather than copy version)
   */
  void setBMatrix(Eigen::Map<siconos::algebra::SiconosSparseMatrix>& S,
                  siconos::algebra::AliasTag);

  /** @brief Set a constant B matrix for the system
   *
   *  Warning : deep copy of the provided object into internal attribute
   *
   *  @param newValue B matrix
   *  @param tag Pass siconos::algebra::copy_t to select this overload (rather than alias
   *
   */
  template <typename T>
  void setBMatrix(T&& newValue, siconos::algebra::CopyTag) {
    static_assert(std::is_same_v<std::decay_t<T>, siconos::algebra::SiconosSparseMatrix>,
                  "Type must be SiconosSparseMatrix");

    if (newValue.rows() != ndof_ || newValue.cols() != ndof_)
      throw std::invalid_argument("Input B matrix has wrong dimensions");

    B_matrix_storage_ =
        std::make_unique<siconos::algebra::SiconosSparseMatrix>(std::forward<T>(newValue));
  }

  /** @return a read-only ref onto the B matrix */
  inline Eigen::Ref<const siconos::algebra::SiconosSparseMatrix> BMatrix() const {
    return use_B_matrix([](Eigen::Ref<const siconos::algebra::SiconosSparseMatrix> M) {
      return Eigen::Ref<const siconos::algebra::SiconosSparseMatrix>(M);
    });
  }

  /*  @return a reference to the B matrix

    pybind11 use only!
  */
  const siconos::algebra::SiconosSparseMatrix& BMatrix_py() const {
    return siconos::algebra::visitStorage<siconos::algebra::AccessMode::OwnedOnly>(
        B_matrix_storage_,
        [](const auto& M) -> const siconos::algebra::SiconosSparseMatrix& { return M; },
        "B_matrix_storage_");
  }

  /** True if B matrix is defined */
  bool hasBMatrix() const {
    return !std::holds_alternative<std::monostate>(B_matrix_storage_);
  }

  inline siconos::algebra::Index stressDimension() const { return dimStress_; }

  /** initialize the internal memory buffers
   *
   *  @param size the size of the siconos::algebra::SiconosMemory. must be >= 0
   */
  void initMemory(siconos::algebra::blocks::size_type size) override;

  inline const siconos::algebra::SiconosMemory& stressMemory() { return stressMemory_; }

  /** push the current values of x, q and r in the stored previous values
   *  xMemory, qMemory, rMemory,
   *  \todo Modify the function swapIn Memory with the new Object Memory
   */
  void swapInMemory() override;

  modeling::Type acceptType(types::FindType& ft) const override { return ft.visit(*this); }
};

}  // namespace siconos::mechanics::fem
#endif  // SOLIDLAGRANGIANTIDS_H
