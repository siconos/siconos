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
#include <vector>

#include "FETypes.hpp"
namespace siconos::mechanics::fem {

class FENode;
class MElement;

using GaussPointsTab = std::vector<std::vector<double>>;

static const GaussPointsTab GaussPointsL2_3 = {
    {1/sqrt(3), 1.0},
    {-1/sqrt(3), 1.0}};

static const GaussPointsTab GaussPointsT3_1 = {
    {0.333333333333333333, 0.333333333333333333, 0.5}};

static const GaussPointsTab GaussPointsT3_2 = {
    {0.66666666666666667, 0.16666666666666667, 0.1666666666666666},
    {0.16666666666666667, 0.66666666666666667, 0.1666666666666666},
    {0.16666666666666667, 0.16666666666666667, 0.1666666666666666}};

static const GaussPointsTab GaussPointsTH4_1 = {{0.25, 0.25, 0.25, 1.0}};

static const GaussPointsTab GaussPointsTH4_2 = {
    {0.13819660112501052, 0.13819660112501052, 0.13819660112501052, 0.04166666666666666},
    {0.13819660112501052, 0.13819660112501052, 0.58541019662496840, 0.04166666666666666},
    {0.13819660112501052, 0.585410196624968405, 0.13819660112501052, 0.04166666666666666},
    {0.58541019662496840, 0.13819660112501052, 0.13819660112501052, 0.04166666666666666}};

static const GaussPointsTab GaussPointsEmpty = {};

/** A Finite Element */
class FElement {
 private:
  /** element number */
  size_t _num = 0;

  /** Element type */
  FiniteElementType _type{FiniteElementType::T3};

  /** Element Family */
  FiniteElementFamily _family{FiniteElementFamily::isoparametric};

  /** number of dof by Element  */
  unsigned _ndof{6};

  /** nodes */
  std::vector<std::shared_ptr<FENode>> _nodes = {};

  /* associated Mesh element */
  std::shared_ptr<MElement> _mElement{nullptr};
  /** Rule of five */
  FElement() = delete;
  FElement(FElement&) = delete;
  FElement& operator=(const FElement&) = delete;
  FElement(FElement&&) = delete;
  FElement& operator=(FElement&&) = delete;

 public:
  /** Constructor
      \param type FE id
      \param ndof dof number
      \param e mesh element
   */
  FElement(FiniteElementType type, unsigned int ndof, std::shared_ptr<MElement> e);

  ~FElement() noexcept = default;

  auto ndof() { return _ndof; }
  auto num() { return _num; }

  auto mElement() { return _mElement; }

  auto family() { return _family; }

  int order();

  int ndofPerNode();

  const GaussPointsTab& GaussPoints(int order);

  std::vector<std::shared_ptr<FENode>>& nodes() { return _nodes; }

  void shapeFunctionIso2D(double ksi, double eta, std::vector<double>& N,
                          std::vector<double>& Nksi, std::vector<double>& Neta);

  void shapeFunctionIso3D(double ksi, double eta, double zeta, std::vector<double>& N,
                          std::vector<double>& Nksi, std::vector<double>& Neta,
                          std::vector<double>& Nzeta);

  void display();
};

}  // namespace siconos::mechanics::fem

#endif
