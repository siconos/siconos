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

#include "MeshUtils.hpp"
#include "Mesh.hpp"
#include <math.h>

#include <fstream>
using namespace std;
#include <stdio.h>
#include <iostream>
#include <sstream>

#include <string>
#include <algorithm>
#include <iterator>
#include <tuple>

#include <stdexcept>

#define DEBUG_STDOUT
#define DEBUG_NOCOLOR
#define DEBUG_MESSAGES
#include "siconos_debug.h"

using namespace siconos::mechanics::fem::native;
template <class Container>
void split(const std::string& str, Container& cont,
           const std::string& delims = " ")
{
  std::size_t current, previous = 0;
  current = str.find_first_of(delims);
  while(current != std::string::npos)
  {
    cont.push_back(str.substr(previous, current - previous));
    previous = current + 1;
    current = str.find_first_of(delims, previous);
  }
  cont.push_back(str.substr(previous, current - previous));
}


Mesh* siconos::mechanics::fem::native::create2dMesh2x1()
{
  MVertex * v1 = new MVertex(1, 0.,0.,0.);
  MVertex * v2 = new MVertex(2, 1.,0.,0.);
  MVertex * v3 = new MVertex(3, 0.,1.,0.);
  MVertex * v4 = new MVertex(4, 1.,1.,0.);

  std::vector<MVertex *> vertices = {v1, v2, v3, v4};


  std::vector<MVertex *> vertices1 = {v1, v2, v3};
  MElement * e1 = new MElement(1, 2, vertices1);

  std::vector<MVertex *> vertices2 = {v2, v4, v3};
  MElement * e2 = new MElement(2, 2, vertices2);

  std::vector<MElement *> elements = {e1, e2};

  return new Mesh(2, vertices, elements);
}

Mesh* siconos::mechanics::fem::native::create2dMeshnxm(int n, int m, double Lx, double Ly)
{

  double lx = Lx/n;
  double ly = Ly/m;


  std::vector<MVertex *> vertices;

  vertices.resize((n+1)*(m+1));

  for(int i =0;  i < n+1; i++)
  {
    for(int j =0;  j < m+1; j++)
    {
      vertices[i+j*(n+1)] = new MVertex(i+j*(n+1), i*lx, j*ly, 0.);
    }
  }
  std::vector<MElement *> elements;
  //elements.resize(2*n*m);
  int element_cnt=0;
  for(int i =0;  i < n; i++)
  {
    for(int j =0;  j < m; j++)
    {
      std::vector<MVertex *> vertices_e_1 = {vertices[i+j*(n+1)], vertices[i+1+(j)*(n+1)], vertices[i+(j+1)*(n+1)]};
      elements.push_back(new MElement(element_cnt++, 2, vertices_e_1));

      std::vector<MVertex *> vertices_e_2 = {vertices[i+1+(j)*(n+1)], vertices[i+1+(j+1)*(n+1)], vertices[i+(j+1)*(n+1)]};
      elements.push_back(new MElement(element_cnt++, 2, vertices_e_2));
    }
  }
  return new Mesh(2, vertices, elements);
}


Mesh* siconos::mechanics::fem::native::createMeshFromGMSH2(std::string gmsh_filename)
{

  ifstream in(gmsh_filename);

  if(!in)
  {
    cout << "Cannot open input file.\n";
    return NULL;
  }

  char str[500];



  float gmsh_version;
  int m =3;
  unsigned int number_of_physical_names=0;
  std::vector< std::tuple<int, string>> physical_entities;
  unsigned int number_of_vertices;
  std::vector<MVertex *> vertices;
  unsigned int number_of_elements;
  std::vector<MElement *> elements;
  double min_z= 1e+64, max_z = -1e+64;
  string toto;

  while(in)
  {
    std::string line;
    std::getline(in, line);
    //if(in) cout << line << endl;
    //std::cout << "line.front()" << line.front() << std::endl;

    // std::vector<std::string> words;
    // split4(line, words);
    // std::copy(words.begin(), words.end(),
    //           std::ostream_iterator<std::string>(std::cout, "\n"));

    if(line.compare("$MeshFormat") == 0)
    {
      std::getline(in, line);
      std::vector<std::string> words;
      split(line, words);
      // std::copy(words.begin(), words.end(),
      //           std::ostream_iterator<std::string>(std::cout, "\n"));

      stringstream token(words[0]);
      token >> gmsh_version;
      std::cout << "gmsh_version : " << gmsh_version << endl;
      if(gmsh_version >= 3.0)
      {
        std::cout << "this simple reader has been written for gmsh v2. Use meshio to translate in gmsh2 format" << std::endl;
        throw std::invalid_argument("Wrong gmsg format");
      }
    }

    if(line.compare("$PhysicalNames") == 0)
    {
      std::getline(in, line);
      stringstream token(line);
      token >> number_of_physical_names;
      std::cout << "number_of_physical_names : " << number_of_physical_names <<  endl;
      physical_entities.resize(number_of_physical_names);
      while(std::getline(in, line))
      {
        if(line.compare("$EndPhysicalNames") == 0) break;
        std::deque<std::string> words;
        split(line, words);
        stringstream t_type(words.front());
        words.pop_front();
        stringstream t_number(words.front());
        words.pop_front();
        std::string name;
        for(std::string s : words)
        {
          name += s + " ";
        }
        name.erase(name.end()-1);
        name.erase(name.end()-1);
        name.erase(name.begin());
        int type, number;
        t_type >> type;
        t_number >> number;
        std::cout << type << " " << number << " name: "<<name<< std::endl;
        physical_entities[number-1] = make_tuple(type, name);
      }

    }

    if(line.compare("$Nodes") == 0)
    {
      std::getline(in, line);
      stringstream token(line);
      token >> number_of_vertices;
      std::cout << "number_of_vertices : " << number_of_vertices <<  endl;
      while(std::getline(in, line))
      {
        if(line.compare("$EndNodes") == 0) break;
        std::vector<std::string> words;
        split(line, words);
        // std::copy(words.begin(), words.end(),
        //           std::ostream_iterator<std::string>(std::cout, "\n"));
        stringstream token(words[0]);
        int vertex_number ;
        token >> vertex_number;
        double x, y, z ;
        stringstream t_x(words[1]);
        t_x >> x;
        stringstream t_y(words[2]);
        t_y >> y;
        stringstream t_z(words[3]);
        t_z >> z;
        max_z = std::max(max_z,z);
        min_z = std::min(min_z,z);
        vertices.push_back(new MVertex(vertex_number, x, y, z));
        //vertices.back()->display();
      }
    }

    if(fabs(min_z-max_z) < 1e-16)
    {
      m =2;
    }
    if(line.compare("$Elements") == 0)
    {
      std::getline(in, line);
      stringstream token(line);
      token >> number_of_elements;
      std::cout << "number_of_elements : " << number_of_elements <<  endl;
      while(std::getline(in, line))
      {
        if(line.compare("$EndElements") == 0) break;
        std::vector<std::string> words;
        split(line, words);
        // std::copy(words.begin(), words.end(),
        //           std::ostream_iterator<std::string>(std::cout, "\n"));
        stringstream token(words[0]);
        int element_number ;
        token >> element_number;
        stringstream t_type(words[1]);
        int element_type ;
        t_type >> element_type;
        stringstream t_nt(words[2]);
        int number_of_tags ;
        t_nt >> number_of_tags;
        std::vector<int> tags ;
        for(int k =3 ; k < 3+number_of_tags; k++)
        {
          int tag;
          stringstream token(words[k]);
          token >> tag;
          tags.push_back(tag);
        }

        std::vector<MVertex *> vertices_e ;
        for(int k = 3+number_of_tags; k < words.size(); k++)
        {
          int node_number;
          stringstream token(words[k]);
          token >> node_number;
          //        std::cout << "node_number" << node_number << std::endl;
          int v;
          for(int k  = 0 ; k <  vertices.size(); k++)
          {
            if(node_number-1+k < vertices.size())
              v  = node_number-1+k;
            else
            {
              v = node_number-1+k -vertices.size();
            }

            if(node_number == vertices[v]->num())
            {
              vertices_e.push_back(vertices[v]);
              break;
            }
          }
        }
        elements.push_back(new MElement(element_number, element_type, vertices_e, tags));
        //elements.back()->display();
      }
    }






  }


  in.close();


  return new Mesh(m, vertices, elements, physical_entities);;
}

void  siconos::mechanics::fem::native::writeMeshforPython(std::shared_ptr<Mesh>  mesh)
{
  FILE * foutput = fopen("mesh.py", "w");
  fprintf(foutput, "coord=[]\n");
  for(MVertex * v : mesh->vertices())
  {
    fprintf(foutput, "coord.append([%e, %e])\n", v->x(), v->y());
  }
  fprintf(foutput, "triangle=[]\n");
  for(MElement * e : mesh->elements())
  {
    fprintf(foutput, "triangle.append([");
    for(MVertex * v : e->vertices())
    {
      fprintf(foutput, "%zu, ", v->num());
    }
    fprintf(foutput, "])\n");
  }
  fclose(foutput);
}
std::string siconos::mechanics::fem::native::prepareWriteDisplacementforPython(std::string basename)
{

  std::string filename = basename + "_displacement.py" ;

  std::cout << "Output displacement for python post-processing in " << filename << std::endl;

  FILE * foutput = fopen(filename.c_str(), "w");
  fprintf(foutput, "import numpy as np\nx=[]\n");
  fprintf(foutput, "y=[]\n");
  fprintf(foutput, "z=[]\n");
  fclose(foutput);
  return filename;
}

void  siconos::mechanics::fem::native::writeDisplacementforPython(std::shared_ptr<Mesh>  mesh,
    std::shared_ptr<FiniteElementModel> femodel,
    std::shared_ptr<SiconosVector> x, std::string filename)
{
  FILE * foutput = fopen(filename.c_str(), "a");
  fprintf(foutput, "x.append(np.array([");

  for(MVertex * v : mesh->vertices())
  {
    std::shared_ptr<FENode> n = femodel->vertexToNode(v);
    double value = 0.0;
    if(n)
    {
      unsigned int idx= (*n->dofIndex())[0];
      value =(*x)(idx);
    }
    fprintf(foutput, "%e,", value) ;

  }
  fprintf(foutput, "]))\n") ;

  fprintf(foutput, "\n");


  fprintf(foutput, "y.append(np.array([");
  for(MVertex * v : mesh->vertices())
  {
    std::shared_ptr<FENode> n = femodel->vertexToNode(v);
    double value = 0.0;
    if(n)
    {
      unsigned int idx= (*n->dofIndex())[1];
      value =(*x)(idx);
    }
    fprintf(foutput, "%e,", value) ;

  }
  fprintf(foutput, "]))\n") ;
  fprintf(foutput, "\n");



  fprintf(foutput, "z.append(np.array([");
  for(MVertex * v : mesh->vertices())
  {
    std::shared_ptr<FENode> n = femodel->vertexToNode(v);
    double value = 0.0;
    if(n)
    {
      if(mesh->dim() > 2)
      {
        unsigned int idx= (*n->dofIndex())[2];
        value =(*x)(idx);
      }
      else
        value=0.0;
    }
    fprintf(foutput, "%e,", value) ;

  }
  fprintf(foutput, "]))\n") ;
  fprintf(foutput, "\n");


  fclose(foutput);
}
