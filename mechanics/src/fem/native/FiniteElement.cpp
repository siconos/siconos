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
#include "Mesh.hpp"  // MElement
#include "SiconosException.hpp"
#include "Tools.hpp"  // enum_to_string

siconos::mechanics::fem::FElement::FElement(FiniteElementType type, unsigned int ndof,
                                            std::shared_ptr<MElement> e)
    : _num(e->num()), _type(type), _ndof(ndof), _mElement(e) {}

int siconos::mechanics::fem::FElement::order() {
  switch (_type) {
    case FiniteElementType::T3:
    case FiniteElementType::TH4:
      return 1;
      break;
    default:
      throw("FElement::order(). element type not recognized");
  }
  return 0;
}

int siconos::mechanics::fem::FElement::ndofPerNode() {
  switch (_type) {
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

double siconos::mechanics::fem::FElement::length() {
  switch (_type) {
    case FiniteElementType::L2:
      return sqrt(pow(_nodes[0]->x()-_nodes[1]->x(),2) + pow(_nodes[0]->y()-_nodes[1]->y(),2) + pow(_nodes[0]->z()-_nodes[1]->z(),2));
    default:
      throw("FElement::norm(). element type not recognized");
  }
}

const siconos::mechanics::fem::GaussPointsTab& siconos::mechanics::fem::FElement::GaussPoints(
    int order) {
  switch (_type) {
    case FiniteElementType::L2:
      return GaussPointsL2_3;
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
  switch (_type) {
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
  switch (_type) {
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
  std::cout << " - FElement - number: " << _num
            << "            - type: " << siconos::tools::enum_to_string(_type)
            << "            - ndof: " << _ndof
            << "            - number of nodes: " << _nodes.size();
  std::cout << std::endl;
  for (std::shared_ptr<FENode> n : _nodes) {
    n->display();
  }
};
