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

#include "FiniteElement.hpp"

#include <iostream>

#include "FENode.hpp"
#include "FETypes.hpp"
#include "Mesh.hpp"  // MElement
#include "SiconosException.hpp"
#include "Tools.hpp"  // enum_to_string

siconos::mechanics::fem::FElement::FElement(FiniteElementType type,
                                            siconos::algebra::Index ndof,
                                            std::shared_ptr<MElement> e)
    : ndof_(ndof), mElement_(e) {}

int siconos::mechanics::fem::FElement::num() const { return mElement_->num(); }

int siconos::mechanics::fem::FElement::order() const {
  FiniteElementType type = mElement_->type();
  switch (type) {
    case FiniteElementType::T3:
    case FiniteElementType::TH4:
      return 1;
      break;
    default:
      throw("FElement::order(). element type not recognized");
  }
  return 0;
}

int siconos::mechanics::fem::FElement::ndofPerNode() const {
  FiniteElementType type = mElement_->type();
  switch (type) {
    case FiniteElementType::T3:
    case FiniteElementType::Q4:
      return 2;
    case FiniteElementType::TH4:
      return 3;
    default:
      throw("FElement::ndorPernode(). element type not recognized");
  }
  return 0;
}

const siconos::mechanics::fem::GaussPointsTab& siconos::mechanics::fem::FElement::GaussPoints(
    int order) {
  FiniteElementType type = mElement_->type();
  switch (type) {
    case FiniteElementType::T3:
      if (order == 1)
        return GaussPointsT3_1;
      else if (order == 2)
        return GaussPointsT3_2;
      break;
    case FiniteElementType::TH4:
      if (order == 1)
        return GaussPointsTH4_1;
      else if (order == 2)
        return GaussPointsTH4_2;
      break;

    default:
      throw("FElement::GaussPoints(). element type not recognized");
  }
  return GaussPointsEmpty;
}

// Note FP : todo, template this class with number of nodes in the element ?
void siconos::mechanics::fem::FElement::shapeFunctionIso2D(double ksi, double eta,
                                                           std::vector<double>& N,
                                                           std::vector<double>& Nksi,
                                                           std::vector<double>& Neta) {
  FiniteElementType type = mElement_->type();
  switch (type) {
    case FiniteElementType::T3: {
      N[0] = 1.0 - ksi - eta;
      N[1] = ksi;
      N[2] = eta;
      Nksi[0] = -1.0;
      Nksi[1] = 1.0;
      Nksi[2] = 0.0;
      Neta[0] = -1.0;
      Neta[1] = 0.0;
      Neta[2] = 1.0;
      break;
    }
    default:
      THROW_EXCEPTION("FElement::shapeFunctionIso2D(). element type not recognized");
  }
}

void siconos::mechanics::fem::FElement::shapeFunctionIso3D(double ksi, double eta, double zeta,
                                                           std::vector<double>& N,
                                                           std::vector<double>& Nksi,
                                                           std::vector<double>& Neta,
                                                           std::vector<double>& Nzeta) {
  FiniteElementType type = mElement_->type();
  switch (type) {
    case FiniteElementType::TH4: {
      N[0] = 1.0 - ksi - eta - zeta;
      N[1] = ksi;
      N[2] = eta;
      N[3] = zeta;

      Nksi[0] = -1.0;
      Nksi[1] = 1.0;
      Nksi[2] = 0.0;
      Nksi[3] = 0.0;

      Neta[0] = -1.0;
      Neta[1] = 0.0;
      Neta[2] = 1.0;
      Neta[3] = 0.0;

      Nzeta[0] = -1.0;
      Nzeta[1] = 0.0;
      Nzeta[2] = 0.0;
      Nzeta[3] = 1.0;

      break;
    }
    default:
      throw("FElement::shapeFunctionIso3D(). element type not recognized");
  }
}

void siconos::mechanics::fem::FElement::display() {
  std::cout << " - FElement - number: " << mElement_->num()
            << "            - type: " << siconos::tools::enum_to_string(mElement_->type())
            << "            - ndof: " << ndof_
            << "            - number of nodes: " << nodes_.size();
  std::cout << std::endl;
  for (std::shared_ptr<FENode> n : nodes_) {
    n->display();
  }
};
