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

#include "Mesh.hpp"

#include <stdio.h>

#include <iostream>
#include <string>
#include <tuple>

#include "Tools.hpp"

void siconos::mechanics::fem::MVertex::display() const {
  std::cout << " - Vertex - number: " << num_ << " ; (x,y,z): " << _x << ", " << _y << ", "
            << _z << "\n";
};

void siconos::mechanics::fem::MElement::display() const {
  std::cout << " - Element - number: " << num_
            << " ; type: " << siconos::tools::enum_to_string(type_) << " ; vertices: ";

  for (auto& v : vertices_) {
    std::cout << " " << v->num();
  }
  std::cout << " - Tags: ";
  for (auto t : tags_) {
    std::cout << " " << t;
  }

  std::cout << "\n";
};

siconos::mechanics::fem::Mesh::Mesh(
    int dim, std::vector<std::shared_ptr<MVertex>> vertices,
    std::vector<std::shared_ptr<MElement>> elements,
    std::vector<std::tuple<int, std::string>> physical_entities)
    : dimension_{dim},
      vertices_{vertices},
      elements_{elements},
      physical_entities_{physical_entities} {
  // Construction of the reverse map : node -> element
  for (auto elem : elements_) {
    for (auto vertex : elem->vertices()) {
      vertex->link_element(elem);
    }
  }
};

void siconos::mechanics::fem::Mesh::display(bool brief) const {
  std::cout << "===== Mesh display ===== \n";
  std::cout << "- dimension : " << dimension_;
  std::cout << "\n- numberOfNodes : " << vertices_.size();
  std::cout << "\n- numberOfElements : " << elements_.size() << "\n";

  int cnt = 0;
  for (const auto& v : vertices_) {
    v->display();
    if (brief and cnt > 10) {
      std::cout << "  ..... \n";
      break;
    }
    cnt++;
  }

  cnt = 0;
  for (auto e : elements_) {
    e->display();
    if (brief and cnt > 10) {
      std::cout << "  ..... \n";
      break;
    }
    cnt++;
  }

  int k = 0;
  std::cout << "Physical entities : (number, tag,  type, name)\n";
  for (auto& physical_entry : physical_entities_) {
    std::cout << k << " " << k + 1 << " " << std::get<0>(physical_entry) << " "
              << std::get<1>(physical_entry) << "\n";
    k++;
  }

  // for(std::vector<MElement *>::iterator it = elements_.begin();
  //     it != elements_.end(); ++it) {
  //   std::cout << *it << std::endl;
  // }
  std::cout << "=========================================================== \n";
}
