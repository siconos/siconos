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

#include "MeshUtils.hpp"

#include <cmath>  // fabs
#include <deque>
#include <fstream>
#include <iostream>
#include <sstream>  // stringstream
#include <string>

#include "FENode.hpp"
#include "FiniteElementModel.hpp"  // FENode
#include "SiconosException.hpp"
#include "SiconosVector.hpp"
#include "RotationQuaternion.hpp"

template <class Container>
void split(const std::string &str, Container &cont, const std::string &delims = " ") {
  std::size_t current, previous = 0;
  current = str.find_first_of(delims);
  while (current != std::string::npos) {
    cont.push_back(str.substr(previous, current - previous));
    previous = current + 1;
    current = str.find_first_of(delims, previous);
  }
  cont.push_back(str.substr(previous, current - previous));
}

std::shared_ptr<siconos::mechanics::fem::Mesh> siconos::mechanics::fem::create2dMesh2x1() {
  auto v1 = std::make_shared<MVertex>(1, 0., 0., 0.);
  auto v2 = std::make_shared<MVertex>(2, 1., 0., 0.);
  auto v3 = std::make_shared<MVertex>(3, 0., 1., 0.);
  auto v4 = std::make_shared<MVertex>(4, 1., 1., 0.);

  std::vector<std::shared_ptr<MVertex>> vertices = {v1, v2, v3, v4};

  std::vector<std::shared_ptr<MVertex>> vertices1 = {v1, v2, v3};
  auto e1 = std::make_shared<MElement>(1, FiniteElementType::T3, vertices1);

  std::vector<std::shared_ptr<MVertex>> vertices2 = {v2, v4, v3};
  auto e2 = std::make_shared<MElement>(2, FiniteElementType::T3, vertices2);

  std::vector<std::shared_ptr<MElement>> elements = {e1, e2};

  return std::make_shared<Mesh>(2, vertices, elements);
}

std::shared_ptr<siconos::mechanics::fem::Mesh> siconos::mechanics::fem::create2dMeshnxm(
    int n, int m, double Lx, double Ly) {
  double lx = Lx / n;
  double ly = Ly / m;

  std::vector<std::shared_ptr<MVertex>> vertices;

  vertices.resize((n + 1) * (m + 1));

  for (int i = 0; i < n + 1; i++) {
    for (int j = 0; j < m + 1; j++) {
      vertices[i + j * (n + 1)] =
          std::make_shared<MVertex>(i + j * (n + 1), i * lx, j * ly, 0.);
    }
  }
  std::vector<std::shared_ptr<MElement>> elements;
  // elements.resize(2*n*m);
  int element_cnt = 0;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      std::vector<std::shared_ptr<MVertex>> vertices_e_1 = {vertices[i + j * (n + 1)],
                                                            vertices[i + 1 + (j) * (n + 1)],
                                                            vertices[i + (j + 1) * (n + 1)]};
      auto e = std::make_shared<MElement>(element_cnt++, FiniteElementType::T3, vertices_e_1);
      elements.push_back(e);

      std::vector<std::shared_ptr<MVertex>> vertices_e_2 = {
          vertices[i + 1 + (j) * (n + 1)], vertices[i + 1 + (j + 1) * (n + 1)],
          vertices[i + (j + 1) * (n + 1)]};

      elements.push_back(
          std::make_shared<MElement>(element_cnt++, FiniteElementType::T3, vertices_e_2));
    }
  }
  return std::make_shared<Mesh>(2, vertices, elements);
}

std::shared_ptr<siconos::mechanics::fem::Mesh> siconos::mechanics::fem::createBeamMesh(std::shared_ptr<MVertex> v0, std::shared_ptr<MVertex> vEnd, int nbBeams, int dim){

  std::vector<std::shared_ptr<MVertex>> vertices;
  vertices.resize(nbBeams+1);
  std::vector<std::shared_ptr<MElement>> elements;

  siconos::mechanics::fem::FiniteElementType FEtype;
  if (dim == 2){
    for (int i=0; i<nbBeams+1;i++){
      vertices[i] = std::make_shared<MVertex>(i+1,i*(vEnd->x() - v0->x())/nbBeams + v0->x(), i*(vEnd->y() - v0->y())/nbBeams + v0->y(), 0.0);
    }
    FEtype = FiniteElementType::B2;
  }
  else if (dim ==3){
    for (int i=0; i<nbBeams+1;i++){
      vertices[i] = std::make_shared<MBeamVertex>(i+1,i*(vEnd->x() - v0->x())/nbBeams + v0->x(), i*(vEnd->y() - v0->y())/nbBeams + v0->y(), i*(vEnd->z() - v0->z())/nbBeams + v0->z(), 0, 0, 0);
    }
    FEtype = FiniteElementType::B3;
  }
  else
    THROW_EXCEPTION("Invalid beam dimension, should be 2 or 3 !");
  std::vector<int> tag = {1};
  std::vector<int> BCtag = {2};
  std::vector<int> Forcetag = {3};
  for (int i = 0; i < nbBeams; i++) {
    std::shared_ptr<MElement> e;
    std::vector<std::shared_ptr<MVertex>> verticesInBeamElement = {vertices[i], vertices[i + 1]};
    if (i == 0)
       e = std::make_shared<MElement>(i+1, FEtype, verticesInBeamElement, BCtag);
    else if (i == nbBeams-1)
       e = std::make_shared<MElement>(i+1, FEtype, verticesInBeamElement, Forcetag);
    else
       e = std::make_shared<MElement>(i+1, FEtype, verticesInBeamElement, tag); // Tag by default
    elements.push_back(e);
  }

  return std::make_shared<Mesh>(dim, vertices, elements);
}

std::shared_ptr<siconos::mechanics::fem::Mesh> siconos::mechanics::fem::createMeshFromGMSH2(
    std::string gmsh_filename) {
  std::ifstream in(gmsh_filename);
  if (!in.is_open()) {
    std::cout << "failed to open " << gmsh_filename << '\n';
    return nullptr;
  }

  // char str[500];
  float gmsh_version;
  int m = 3;
  unsigned int number_of_physical_names = 0;
  std::vector<std::tuple<int, std::string>> physical_entities;
  unsigned int number_of_vertices;
  std::vector<std::shared_ptr<MVertex>> vertices;
  unsigned int number_of_elements;
  std::vector<std::shared_ptr<MElement>> elements;
  double min_z = 1e+64, max_z = -1e+64;
  std::string toto;

  while (in) {
    std::string line;
    std::getline(in, line);
    if (line.compare("$MeshFormat") == 0) {
      std::getline(in, line);
      std::vector<std::string> words;
      split(line, words);
      std::stringstream token(words[0]);
      token >> gmsh_version;
      std::cout << "gmsh_version : " << gmsh_version << "\n";
      if (gmsh_version >= 3.0) {
        std::cout << "this simple reader has been written for gmsh v2. Use meshio to "
                     "translate in gmsh2 format\n";
        THROW_EXCEPTION("Wrong gmsh format");
      }
    }

    if (line.compare("$PhysicalNames") == 0) {
      std::getline(in, line);
      std::stringstream token(line);
      token >> number_of_physical_names;
      std::cout << "number_of_physical_names : " << number_of_physical_names << "\n";
      physical_entities.resize(number_of_physical_names);
      while (std::getline(in, line)) {
        if (line.compare("$EndPhysicalNames") == 0) break;
        std::deque<std::string> words;
        split(line, words);
        std::stringstream t_type(words.front());
        words.pop_front();
        std::stringstream t_number(words.front());
        words.pop_front();
        std::string name;
        for (std::string s : words) {
          name += s + " ";
        }
        name.erase(name.end() - 1);
        name.erase(name.end() - 1);
        name.erase(name.begin());
        int type, number;
        t_type >> type;
        t_number >> number;
        std::cout << type << " " << number << " name: " << name << std::endl;
        physical_entities[number - 1] = make_tuple(type, name);
      }
    }

    if (line.compare("$Nodes") == 0) {
      std::getline(in, line);
      std::stringstream token(line);
      token >> number_of_vertices;
      std::cout << "number_of_vertices : " << number_of_vertices << "\n";
      while (std::getline(in, line)) {
        if (line.compare("$EndNodes") == 0) break;
        std::vector<std::string> words;
        split(line, words);
        // std::copy(words.begin(), words.end(),
        //           std::ostream_iterator<std::string>(std::cout, "\n"));
        std::stringstream token(words[0]);
        int vertex_number;
        token >> vertex_number;
        double x, y, z;
        std::stringstream t_x(words[1]);
        t_x >> x;
        std::stringstream t_y(words[2]);
        t_y >> y;
        std::stringstream t_z(words[3]);
        t_z >> z;
        max_z = std::max(max_z, z);
        min_z = std::min(min_z, z);
        vertices.push_back(std::make_shared<MVertex>(vertex_number, x, y, z));
        // vertices.back()->display();
      }
    }

    if (std::fabs(min_z - max_z) < 1e-16) {
      m = 2;
    }
    if (line.compare("$Elements") == 0) {
      std::getline(in, line);
      std::stringstream token(line);
      token >> number_of_elements;
      std::cout << "number_of_elements : " << number_of_elements << "\n";
      while (std::getline(in, line)) {
        if (line.compare("$EndElements") == 0) break;
        std::vector<std::string> words;
        split(line, words);
        // std::copy(words.begin(), words.end(),
        //           std::ostream_iterator<std::string>(std::cout, "\n"));
        std::stringstream token(words[0]);
        int element_number;
        token >> element_number;
        std::stringstream t_type(words[1]);
        int tempvalue;
        t_type >> tempvalue;
        FiniteElementType element_type =
            (FiniteElementType)tempvalue;  // This has to be checked carefully ...

        std::stringstream t_nt(words[2]);
        decltype(words.size()) number_of_tags;
        t_nt >> number_of_tags;
        std::vector<int> tags;
        for (decltype(words.size()) k = 3; k < 3 + number_of_tags; k++) {
          int tag;
          std::stringstream token(words[k]);
          token >> tag;
          tags.push_back(tag);
        }

        std::vector<std::shared_ptr<MVertex>> vertices_e;
        for (auto k = 3 + number_of_tags; k < words.size(); k++) {
          decltype(vertices.size()) node_number;
          std::stringstream token(words[k]);
          token >> node_number;
          //        std::cout << "node_number" << node_number << std::endl;
          int v;
          for (decltype(vertices.size()) k = 0; k < vertices.size(); k++) {
            if (node_number - 1 + k < vertices.size())
              v = node_number - 1 + k;
            else {
              v = node_number - 1 + k - vertices.size();
            }

            if (node_number == vertices[v]->num()) {
              vertices_e.push_back(vertices[v]);
              break;
            }
          }
        }
        elements.push_back(
            std::make_shared<MElement>(element_number, element_type, vertices_e, tags));
        // elements.back()->display();
      }
    }
  }

  in.close();

  return std::make_shared<Mesh>(m, vertices, elements, physical_entities);
}

void siconos::mechanics::fem::writeMeshforPython(std::shared_ptr<Mesh> mesh) {
  FILE *foutput = fopen("mesh.py", "w");
  fprintf(foutput, "coord=[]\n");
  for (auto v : mesh->vertices()) {
    fprintf(foutput, "coord.append([%e, %e])\n", v->x(), v->y());
  }
  fprintf(foutput, "triangle=[]\n");
  for (auto e : mesh->elements()) {
    fprintf(foutput, "triangle.append([");
    for (auto v : e->vertices()) {
      fprintf(foutput, "%zu, ", v->num());
    }
    fprintf(foutput, "])\n");
  }
  fclose(foutput);
}
std::string siconos::mechanics::fem::prepareWriteDisplacementforPython(std::string basename) {
  std::string filename = basename + "_displacement.py";

  std::cout << "Output displacement for python post-processing in " << filename << std::endl;

  FILE *foutput = fopen(filename.c_str(), "w");
  fprintf(foutput, "import numpy as np\nx=[]\n");
  fprintf(foutput, "y=[]\n");
  fprintf(foutput, "z=[]\n");
  fclose(foutput);
  return filename;
}
std::string siconos::mechanics::fem::prepareWriteTensorforPython(std::string basename, std::string tensorName) {
  std::string filename = basename + "_tensor.py";

  std::cout << "Output Tensor for python post-processing in " << filename << std::endl;

  FILE *foutput = fopen(filename.c_str(), "w");
  fprintf(foutput, "import numpy as np\n");

  std::string strxx = tensorName + "_xx=[]\n";
  fprintf(foutput, strxx.c_str());

  std::string stryy = tensorName + "_yy=[]\n";
  fprintf(foutput, stryy.c_str());

  std::string strxy = tensorName + "_xy=[]\n";
  fprintf(foutput, strxy.c_str());

  fclose(foutput);
  return filename;
}

std::string siconos::mechanics::fem::prepareWriteBeamTensorforPython(std::string basename, std::string tensorName) {
  std::string filename = basename + "_tensor.py";

  std::cout << "Output Tensor for python post-processing in " << filename << std::endl;

  FILE *foutput = fopen(filename.c_str(), "w");
  fprintf(foutput, "import numpy as np\n");

  std::string strTC = tensorName + "_tension=[]\n";
  fprintf(foutput, strTC.c_str());

  std::string strBending1 = tensorName + "_bending1=[]\n";
  fprintf(foutput, strBending1.c_str());

  std::string strBending2 = tensorName + "_bending2=[]\n";
  fprintf(foutput, strBending2.c_str());

  fclose(foutput);
  return filename;
}

void siconos::mechanics::fem::writeDisplacementforPython(
    std::shared_ptr<Mesh> mesh, std::shared_ptr<FiniteElementModel> femodel,
    std::shared_ptr<siconos::algebra::SiconosVector> x, std::string filename) {
  FILE *foutput = fopen(filename.c_str(), "a");
  fprintf(foutput, "x.append(np.array([");

  for (auto v : mesh->vertices()) {
    auto n = femodel->vertexToNode(v);
    double value = 0.0;
    if (n) {
      unsigned int idx = (*n->dofIndex())[0];
      value = (*x)(idx);
    }
    fprintf(foutput, "%e,", value);
  }
  fprintf(foutput, "]))\n");

  fprintf(foutput, "\n");

  fprintf(foutput, "y.append(np.array([");
  for (auto v : mesh->vertices()) {
    auto n = femodel->vertexToNode(v);
    double value = 0.0;
    if (n) {
      auto idx = (*n->dofIndex())[1];
      value = (*x)(idx);
    }
    fprintf(foutput, "%e,", value);
  }
  fprintf(foutput, "]))\n");
  fprintf(foutput, "\n");

  fprintf(foutput, "z.append(np.array([");
  for (auto v : mesh->vertices()) {
    std::shared_ptr<FENode> n = femodel->vertexToNode(v);
    double value = 0.0;
    if (n) {
      if (mesh->dim() > 2) {
        auto idx = (*n->dofIndex())[2];
        value = (*x)(idx);
      } else
        value = 0.0;
    }
    fprintf(foutput, "%e,", value);
  }
  fprintf(foutput, "]))\n");
  fprintf(foutput, "\n");

  fclose(foutput);
}

void siconos::mechanics::fem::writeTensorforPython(
     std::shared_ptr<FiniteElementModel> femodel,
    std::shared_ptr<siconos::algebra::SiconosVector> x, std::string filename, std::string tensorName) {
  FILE *foutput = fopen(filename.c_str(), "a");

  std::string stringxx = tensorName + "_xx.append(np.array([";
  x->display();
  fprintf(foutput, stringxx.c_str());

  unsigned int elem_cnt = 0;
  double value;
  for(std::shared_ptr<FElement> fe : femodel->elements())
  {
    value = (*x)(elem_cnt*3);
    fprintf(foutput, "%e,", value);
    elem_cnt++;
  }
  std::cout << " elem_cnt:" << elem_cnt << std::endl;

  fprintf(foutput, "]))\n");
  fprintf(foutput, "\n");

  std::string stringyy = tensorName + "_yy.append(np.array([";
  fprintf(foutput, stringyy.c_str());
  elem_cnt = 0;
  for(std::shared_ptr<FElement> fe : femodel->elements())
  {
    value = (*x)(elem_cnt*3+1);
    fprintf(foutput, "%e,", value);
    elem_cnt++;
  }

  fprintf(foutput, "]))\n");
  fprintf(foutput, "\n");

  std::string stringxy = tensorName + "_xy.append(np.array([";
  fprintf(foutput, stringxy.c_str());
  elem_cnt = 0;
  for(std::shared_ptr<FElement> fe : femodel->elements())
  {
    value = (*x)(elem_cnt*3+2);
    fprintf(foutput, "%e,", value);
    elem_cnt++;
  }

  fprintf(foutput, "]))\n");
  fprintf(foutput, "\n");
  fclose(foutput);
}



void siconos::mechanics::fem::prepareWriteBeamPositionforSOFA(std::string filename) {

  std::cout << "Output displacement for SOFA post-processing in " << filename << std::endl;

  FILE *foutput = fopen(filename.c_str(), "w");
  fclose(foutput);
}

void siconos::mechanics::fem::prepareWriteBlockPositionforSOFA(std::string filename) {

  std::cout << "Output displacement for SOFA post-processing in " << filename << std::endl;

  FILE *foutput = fopen(filename.c_str(), "w");
  fclose(foutput);
}

void  siconos::mechanics::fem::writeBeamPositionforSOFA(std::shared_ptr<Mesh>  mesh,
                                                               std::shared_ptr<FiniteElementModel> femodel,
                                                               std::shared_ptr<siconos::algebra::SiconosVector> x, std::string filename, double t)
{
  FILE * foutput = fopen(filename.c_str(), "a");
  fprintf(foutput, "T= %f\n", t);
  auto vertices =  mesh->vertices();
  fprintf(foutput, "X=");
  for(auto v : mesh->vertices())
  {
    std::shared_ptr<FENode> n = femodel->vertexToNode(v);
    double value = 0.0, valueqx, valueqy, valueqz;
    if(n)
    {
      unsigned int idx, idy, idtheta,idz, idqx, idqy, idqz;
      idx = (*n->dofIndex())[0];
      idy = (*n->dofIndex())[1];

      value =(*x)(idx);
      fprintf(foutput, " %e", value+v->x()) ;

      value =(*x)(idy);
      fprintf(foutput, " %e", value+v->y()) ;

      if ((*n->dofIndex()).size() == 3){
        // 2D case
        fprintf(foutput, " %e", 0.0) ;

        idtheta = (*n->dofIndex())[2];
        // idqx= (*n->dofIndex())[3];
        valueqx =(*x)(idtheta);
        // idqy= (*n->dofIndex())[4];
        // valueqy =(*x)(idqy);
        // idqz= (*n->dofIndex())[5];
        // valueqz =(*x)(idqz);

        double quat[4];

        double c3 = cos( (valueqx+v->z()) / 2 );
        double c1 = cos( 0 );
        double c2 = cos( 0 );

        double s3 = sin( (valueqx+v->z()) / 2 );
        double s1 = sin( 0 );
        double s2 = sin( 0 );

        quat[0] = s1 * c2 * c3 - c1 * s2 * s3;
        quat[1] = c1 * s2 * c3 + s1 * c2 * s3;
        quat[2] = c1 * c2 * s3 - s1 * s2 * c3;
        quat[3] = c1 * c2 * c3 + s1 * s2 * s3;

        fprintf(foutput, " %e", quat[0]) ;
        fprintf(foutput, " %e", quat[1]) ;
        fprintf(foutput, " %e", quat[2]) ;
        fprintf(foutput, " %e", quat[3]) ;
      }
      else
      {
        idz = (*n->dofIndex())[2];
        value =(*x)(idz);
        fprintf(foutput, " %e", value+v->z()) ;

        idqx = (*n->dofIndex())[3];
        valueqx =(*x)(idqx);
        idqy = (*n->dofIndex())[4];
        valueqy =(*x)(idqy);
        idqz = (*n->dofIndex())[5];
        valueqz =(*x)(idqz);

        auto angle = std::make_shared<siconos::algebra::SiconosVector>(3);
        angle->setValue(0,valueqx);
        angle->setValue(1,valueqy);
        angle->setValue(2,valueqz);
        auto quat = std::make_shared<siconos::algebra::SiconosVector>(7);

        siconos::geometry::quaternionFromRotationVector(angle,quat);

        fprintf(foutput, " %e", (*quat)(4)) ;
        fprintf(foutput, " %e", (*quat)(5)) ;
        fprintf(foutput, " %e", (*quat)(6)) ;
        fprintf(foutput, " %e", (*quat)(3)) ;
      }

    }
  }
  fprintf(foutput, "\n");
  fclose(foutput);
}

void  siconos::mechanics::fem::writeBlockPositionforSOFA(std::shared_ptr<siconos::algebra::SiconosVector> x, std::string filename, double t)
{
  FILE * foutput = fopen(filename.c_str(), "a");
  fprintf(foutput, "T= %f\n", t);
  fprintf(foutput, "X=");
  fprintf(foutput, " %e", (*x)(0)) ;
  fprintf(foutput, " %e", (*x)(1)) ;
  // 2D case
  fprintf(foutput, " %e", 0.0) ;

  fprintf(foutput, "\n");
  fclose(foutput);
}

