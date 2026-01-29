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
#ifndef FEHPP
#define FEHPP

#include <cmath>
#include <memory>
#include <span>
#include <unordered_set>
#include <vector>

#include "FETypes.hpp"
#include "SiconosMatrix.hpp"
namespace siconos::mechanics::fem {

class FENode;
class MElement;
class MVertex;

using GaussPointsTab = std::vector<std::vector<double>>;

// B2 / B3, order 3
static const GaussPointsTab GaussPointsB2_3 = {{1 / sqrt(3), 1.0}, {-1 / sqrt(3), 1.0}};

// T3, order 1
static const GaussPointsTab GaussPointsT3_1 = {
    {0.333333333333333333, 0.333333333333333333, 0.5}};

// T3, order 2
static const GaussPointsTab GaussPointsT3_2 = {
    {0.66666666666666667, 0.16666666666666667, 0.1666666666666666},
    {0.16666666666666667, 0.66666666666666667, 0.1666666666666666},
    {0.16666666666666667, 0.16666666666666667, 0.66666666666666667}};

// TH4 order 1
static const GaussPointsTab GaussPointsTH4_1 = {{0.25, 0.25, 0.25, 1.0}};

// TH4 order 2
static const GaussPointsTab GaussPointsTH4_2 = {
    {0.13819660112501052, 0.13819660112501052, 0.13819660112501052, 0.04166666666666666},
    {0.13819660112501052, 0.13819660112501052, 0.58541019662496840, 0.04166666666666666},
    {0.13819660112501052, 0.585410196624968405, 0.13819660112501052, 0.04166666666666666},
    {0.58541019662496840, 0.13819660112501052, 0.13819660112501052, 0.04166666666666666}};

static const GaussPointsTab GaussPointsEmpty = {};

/** A Finite Element */
class FElement {
 protected:
  /** element number */
  size_t num_ = 0;

  /** Type */
  FiniteElementType type_{FiniteElementType::Unknown};

  /** Element Family */
  FiniteElementFamily family_{FiniteElementFamily::isoparametric};

  /** number of dof by Element  */
  siconos::algebra::Index ndof_{6};

  /** stress dimension */
  siconos::algebra::Index dimStress_{0};

  /** nodes */
  std::vector<std::shared_ptr<FENode>> nodes_ = {};

  /** Global To Local Matrix */
  siconos::algebra::SiconosDenseMatrix Te_{};  // 6x6 (dim=2) or 12x12(dim=3)

  /** Transpose of Te */
  siconos::algebra::SiconosDenseMatrix Te_transpose_{};  // 6x6 (dim=2) or 12x12(dim=3)

  /* associated Mesh element */
  std::shared_ptr<MElement> mElement_{nullptr};

  /** number of dof per node */
  int ndofPernode_{0};

  /** Element length */
  double length_{0};

  /** Rule of five */
  FElement() = delete;
  FElement(FElement&) = delete;
  FElement& operator=(const FElement&) = delete;
  FElement(FElement&&) = delete;
  FElement& operator=(FElement&&) = delete;

  /** Build Te, Te^t and compute length - Only for B2 and B3 */
  void initialize_beam_element();

  /** @return true if the current element is a valid beam element */
  inline bool is_valid_beam_element(FiniteElementType value) {
    static const std::unordered_set<FiniteElementType> valid_values = {FiniteElementType::B2,
                                                                       FiniteElementType::B3};
    return valid_values.count(value) > 0;
  }

 public:
  /** Constructor
   *  @param elem mesh element used as source of this finite element
   *  @param nodes list of FENode of the current element (built from vertices of the
   *    mesh_element)
   */
  FElement(std::shared_ptr<MElement> e, const std::vector<std::shared_ptr<FENode>>& nodes);

  ~FElement() noexcept = default;

  siconos::algebra::Index ndof() const { return ndof_; }

  siconos::algebra::Index dimStress() const { return dimStress_; }

  auto num() const { return num_; };

  auto type() const { return (type_); }

  void add_node(std::shared_ptr<MVertex> vertex);

  std::shared_ptr<MElement> mElement() { return mElement_; }

  FiniteElementFamily family() const { return family_; }

  int order() const;

  int ndofPerNode() const { return number_of_dof_per_node(type_); };

  std::span<const std::vector<double>> GaussPoints(int order) const;

  //  std::vector<std::shared_ptr<FENode>>& nodes() { return nodes_; }
  std::span<const std::shared_ptr<FENode>> nodes() const { return nodes_; }

  void shapeFunctionIso2D(double ksi, double eta, std::vector<double>& N,
                          std::vector<double>& Nksi, std::vector<double>& Neta) const;

  void shapeFunctionIso3D(double ksi, double eta, double zeta, std::vector<double>& N,
                          std::vector<double>& Nksi, std::vector<double>& Neta,
                          std::vector<double>& Nzeta) const;

  /** @return element length
   *
   *  Only for B2 and B3
   */
  double length() const { return length_; };

  /** @Returns a read-only view onto Te */
  Eigen::Ref<const siconos::algebra::SiconosDenseMatrix> TeMatrix() const { return Te_; }

  /** Returns a read-only view onto transpose(Te) */
  Eigen::Ref<const siconos::algebra::SiconosDenseMatrix> TeMatrix_transpose() const {
    return Te_transpose_;
  }

  void display();
};

}  // namespace siconos::mechanics::fem

#endif
