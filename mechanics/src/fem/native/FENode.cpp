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

#include "FENode.hpp"

#include <iostream>

#include "Mesh.hpp"

double siconos::mechanics::fem::FENode::x() { return _mVertex->x(); }

double siconos::mechanics::fem::FENode::y() { return _mVertex->y(); }

double siconos::mechanics::fem::FENode::z() { return _mVertex->z(); }

void siconos::mechanics::fem::FENode::display() {
  std::cout << "     - Fe Node - number: " << _num
            << "               - ndof:" << _dofIndex->size()
            << "               - dofIndex: " << _dofIndex->front() << ":" << _dofIndex->back() << "\n";
};
