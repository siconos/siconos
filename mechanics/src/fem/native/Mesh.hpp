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

/*! \file Mesh.hpp

 */
#ifndef MESH_H
#define MESH_H

#include <vector>
#include <iostream>
#include <string>



namespace siconos::mechanics::fem::native
{
/** a Mesh container
 */

class MElement;
// a Mesh vertex
struct MVertex
{
  /* Vertex Number */
  size_t _num;

  /* Vextex Coordinate */
  double _x;
  double _y;
  double _z;

  /** elements to which the node belongs **/
  std::vector<MElement *> _elements;


  /* Constructor from data */
  MVertex(size_t num, double x, double y, double z):_num(num), _x(x), _y(y), _z(z) {};

  double x()
  {
    return _x;
  };
  double y()
  {
    return _y;
  };
  double z()
  {
    return _z;
  };


  size_t num()
  {
    return _num;
  }


  std::vector<MElement *> & elements()
  {
    return _elements;
  };

  void display()
  {
    std::cout << " - Vertex - number: " << _num
              << " ; (x,y,z): "
              << _x <<", "
              << _y <<", "
              << _z  ;
    std::cout << std::endl;
  };
};


// a mesh element
class MElement
{
protected :

  /** Element number */
  size_t _num;

  /** type (following gmsh convention) */
  int _type;

  /** vertices **/
  std::vector<MVertex *> _vertices;

  /** tags **/
  std::vector<int> _tags;



  /** default constructor */
  MElement() {};

public:
  MElement(size_t num, int type, std::vector<MVertex *> vertices):
    _num(num),_type(type),_vertices(vertices)
  {
    _tags.push_back(0);
  };

  MElement(size_t num, int type, std::vector<MVertex *> vertices, std::vector<int> tags):
    _num(num),_type(type),_vertices(vertices),_tags(tags) {};

  int type()
  {
    return _type;
  }
  size_t num()
  {
    return _num;
  }
  std::vector<MVertex *>& vertices()
  {
    return _vertices;
  }
  int tags(int n)
  {
    return _tags[n];
  }
  void display()
  {
    std::cout << " - Element - number: " << _num
              << " ; type: " << _type
              << " ; vertices: ";

    for(MVertex * v :_vertices)
    {
      std::cout << " " << v->_num ;
    }
    std::cout << " - Tags: ";
    for(int  t :_tags)
    {
      std::cout << " " << t ;
    }

    std::cout << std::endl;
  };
};


class Mesh
{

protected:
  /* serialization hooks */
  //ACCEPT_SERIALIZATION(Mesh);

  /** Space dimension */
  int _dim;

  /** number of nodes */
  int _numberOfVertices;

  /** number of elements */
  int _numberOfElements;

  /** vertices */
  std::vector<MVertex *> _vertices;

  /** elements */
  std::vector<MElement *> _elements;

  /** Physical entities
   * This vector enables to link the tags to Physical entities
   */
  std::vector< std::tuple<int, std::string>> _physical_entities;


  /** default constructor */
  Mesh() {};

public:
  Mesh(int dim, int numberOfNodes, int numberOfElements);

  /** constructor
   *  \param dim dimension
   */
  Mesh(int dim,
       std::vector<MVertex *> vertices,
       std::vector<MElement *> elements);

  /** constructor
   *  \param dim dimension
   */
  Mesh(int dim,
       std::vector<MVertex *> vertices,
       std::vector<MElement *> elements,
       std::vector< std::tuple<int, std::string>> physical_entities);

  /** destructor */
  ~Mesh() {};

  int dim()
  {
    return _dim;
  };

  std::vector<MVertex *> & vertices()
  {
    return _vertices;
  }

  std::vector<MElement *> & elements()
  {
    return _elements;
  }

  std::vector< std::tuple<int, std::string>> physical_entities()
  {
    return _physical_entities;
  }
  /** print the data of the Mesh
   */
  void display(bool brief = true) const;

  //
  //ACCEPT_STD_VISITORS();

};


} // namespace siconos::mechanics::fem::native
#endif // MESH_H
