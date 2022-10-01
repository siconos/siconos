/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2020 INRIA.
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

#define DEBUG_STDOUT
#define DEBUG_NOCOLOR
#define DEBUG_MESSAGES
#include "siconos_debug.h"

#include <stdio.h>
#include <iostream>
#include <string>
// --- Constructor from a list of attributes
siconos::mechanics::fem::native::Mesh::Mesh(int dim, int numberOfVertices, int numberOfElements):
  _dim(dim),_numberOfVertices(numberOfVertices),_numberOfElements(numberOfElements)
{
  DEBUG_BEGIN("siconos::mechanics::fem:native::Mesh::Mesh(int dim, int numberOfNodes, int numberOfElements )\n");
  DEBUG_END("siconos::mechanics::fem:native::Mesh::Mesh(int dim, int numberOfNodes, int numberOfElements) \n");
};

siconos::mechanics::fem::native::Mesh::Mesh(int dim,
    std::vector<MVertex *> vertices,
    std::vector<MElement *> elements): _dim(dim), _vertices(vertices), _elements(elements)
{
  _numberOfElements =  _elements.size();
  _numberOfVertices =  _vertices.size();

  // Construction of the reverse map : node -> element
  for(MElement * e : _elements)
  {
    for(MVertex * v : e->vertices())
    {
      v->elements().push_back(e);
    }
  }

  _physical_entities= {};

};

siconos::mechanics::fem::native::Mesh::Mesh(int dim,
    std::vector<MVertex *> vertices,
    std::vector<MElement *> elements,
    std::vector< std::tuple<int, std::string>> physical_entities): Mesh(dim, vertices, elements)
{
  _physical_entities=physical_entities;
};

void siconos::mechanics::fem::native::Mesh::display(bool brief) const
{
  std::cout << "===== Mesh display ===== " <<std::endl;
  std::cout << "- dimension : " << _dim << std::endl;
  std::cout << "- numberOfNodes : " << _numberOfVertices << std::endl;
  std::cout << "- numberOfElements : " << _numberOfElements << std::endl;


  int cnt =0;
  for(MVertex * v : _vertices)
  {
    v->display();
    if(brief and cnt >10)
    {
      std::cout << "  ..... " << std::endl;
      break;
    }
    cnt++;
  }

  cnt=0;
  for(MElement * e : _elements)
  {
    e->display();
    if(brief and cnt >10)
    {
      std::cout << "  ..... " << std::endl;
      break;
    }
    cnt++;
  }

  int k=0;
  std::cout << "Physical entities : (number, tag,  type, name)" << std::endl;
  for(std::tuple<int, std::string> physical_entry :  _physical_entities)
  {
    std::cout << k  << " "  <<  k+1 << " " << std::get<0>(physical_entry) << " " <<  std::get<1>(physical_entry) <<    std::endl;
    k++;
  }


  // for(std::vector<MElement *>::iterator it = _elements.begin();
  //     it != _elements.end(); ++it) {
  //   std::cout << *it << std::endl;
  // }
  std::cout << "=========================================================== " <<std::endl;

}
