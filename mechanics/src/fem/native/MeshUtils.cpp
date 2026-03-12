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
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>  // stringstream
#include <string>

#include "FENode.hpp"
#include "FETypes.hpp"
#include "FemTools.hpp"
#include "FiniteElementModel.hpp"
#include "Mesh.hpp"  // For MeshVertex, MeshElement ...
#include "RotationQuaternion.hpp"
#include "SiconosException.hpp"
#include "SiconosVector.hpp"

template <class Container>
void split(const std::string& str, Container& cont, const std::string& delims = " ") {
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
  auto v1 = std::make_shared<MeshVertex>(1, 0., 0., 0.);
  auto v2 = std::make_shared<MeshVertex>(2, 1., 0., 0.);
  auto v3 = std::make_shared<MeshVertex>(3, 0., 1., 0.);
  auto v4 = std::make_shared<MeshVertex>(4, 1., 1., 0.);

  std::vector<std::shared_ptr<MeshVertex>> vertices = {v1, v2, v3, v4};

  std::vector<std::shared_ptr<MeshVertex>> vertices1 = {v1, v2, v3};
  auto e1 = std::make_shared<MeshElement>(1, FiniteElementType::T3, vertices1);

  std::vector<std::shared_ptr<MeshVertex>> vertices2 = {v2, v4, v3};
  auto e2 = std::make_shared<MeshElement>(2, FiniteElementType::T3, vertices2);

  std::vector<std::shared_ptr<MeshElement>> elements = {e1, e2};

  return std::make_shared<Mesh>(2, vertices, elements);
}

std::shared_ptr<siconos::mechanics::fem::Mesh> siconos::mechanics::fem::create2dMeshnxm(
    size_t n, size_t m, double Lx, double Ly) {
  double lx = Lx / static_cast<double>(n);
  double ly = Ly / static_cast<double>(m);

  std::vector<std::shared_ptr<MeshVertex>> vertices;

  vertices.resize((n + 1) * (m + 1));

  for (size_t i = 0; i < n + 1; i++) {
    for (size_t j = 0; j < m + 1; j++) {
      vertices[i + j * (n + 1)] = std::make_shared<MeshVertex>(
          i + j * (n + 1), static_cast<double>(i) * lx, static_cast<double>(j) * ly, 0.);
    }
  }
  std::vector<std::shared_ptr<MeshElement>> elements;
  // elements.resize(2*n*m);
  int element_cnt = 0;
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < m; j++) {
      std::vector<std::shared_ptr<MeshVertex>> vertices_e_1 = {
          vertices[i + j * (n + 1)], vertices[i + 1 + (j) * (n + 1)],
          vertices[i + (j + 1) * (n + 1)]};
      auto e =
          std::make_shared<MeshElement>(element_cnt++, FiniteElementType::T3, vertices_e_1);
      elements.push_back(e);

      std::vector<std::shared_ptr<MeshVertex>> vertices_e_2 = {
          vertices[i + 1 + (j) * (n + 1)], vertices[i + 1 + (j + 1) * (n + 1)],
          vertices[i + (j + 1) * (n + 1)]};

      elements.push_back(
          std::make_shared<MeshElement>(element_cnt++, FiniteElementType::T3, vertices_e_2));
    }
  }
  return std::make_shared<Mesh>(2, vertices, elements);
}

std::shared_ptr<siconos::mechanics::fem::Mesh> siconos::mechanics::fem::createBeamMesh(
    const siconos::algebra::SiconosVector3& coords_start,
    const siconos::algebra::SiconosVector3& coords_end, size_t nb_elements, int dim) {
  // 1. Computes the list of vertices
  siconos::algebra::SiconosVector3 step = (coords_end - coords_start) / nb_elements;

  siconos::mechanics::fem::FiniteElementType FEtype;
  if (dim == 2) {
    FEtype = FiniteElementType::B2;
    step(2) = 0;
  } else if (dim == 3)
    FEtype = FiniteElementType::B3;
  else
    THROW_EXCEPTION("Invalid beam dimension, should be 2 or 3!");

  std::vector<std::shared_ptr<MeshVertex>> vertices(nb_elements + 1);

  for (size_t i = 0; i < nb_elements + 1; i++) {
    vertices[i] = std::make_shared<MeshBeamVertex>(
        i, coords_start[0] + step(0) * static_cast<double>(i),
        coords_start[1] + static_cast<double>(i) * step(1),
        coords_start[2] + static_cast<double>(i) * step(2), 0, 0, 0);
  }

  // 2. Computes the list of mesh elements
  // We assume:
  // - boundary condition tag for the first vertex/element
  // - applied force tag for the last element
  // - default tag for the others
  std::vector<std::shared_ptr<MeshElement>> elements;

  // First element
  std::vector<std::shared_ptr<MeshVertex>> vertices_in_elem = {vertices[0], vertices[1]};

  std::vector<std::shared_ptr<MeshVertex>> vertices_start = {vertices[0], vertices[1]};
  std::vector<int> tags = {
      default_tags.at(siconos::mechanics::fem::MeshTags::boundary_conditions)};
  elements.push_back(std::make_shared<MeshElement>(0, FEtype, vertices_start, tags));
  for (size_t i = 1; i < nb_elements - 1; i++) {
    tags = {default_tags.at(siconos::mechanics::fem::MeshTags::bulk_material)};
    std::vector<std::shared_ptr<MeshVertex>> vertices_in_elem = {vertices[i], vertices[i + 1]};
    elements.push_back(std::make_shared<MeshElement>(i, FEtype, vertices_in_elem, tags));
  }
  tags = {default_tags.at(siconos::mechanics::fem::MeshTags::applied_forces)};
  std::vector<std::shared_ptr<MeshVertex>> vertices_in_last_elem = {vertices[nb_elements - 1],
                                                                    vertices[nb_elements]};
  elements.push_back(
      std::make_shared<MeshElement>(nb_elements - 1, FEtype, vertices_in_last_elem, tags));
  elements.back()->display();
  return std::make_shared<Mesh>(dim, vertices, elements);
}

std::shared_ptr<siconos::mechanics::fem::Mesh> siconos::mechanics::fem::createMeshFromGMSH2(
    std::string gmsh_filename) {
  std::ifstream in(gmsh_filename);
  if (!in.is_open()) {
    throw std::runtime_error("Cannot open file: " + gmsh_filename);
  }

  // char str[500];
  float gmsh_version;
  int m = 3;
  unsigned int number_of_physical_names = 0;
  std::vector<std::tuple<int, std::string>> physical_entities;
  unsigned int number_of_vertices;
  std::vector<std::shared_ptr<MeshVertex>> vertices;
  unsigned int number_of_elements;
  std::vector<std::shared_ptr<MeshElement>> elements;
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
        int type;
        size_t number;
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
        vertices.push_back(std::make_shared<MeshVertex>(vertex_number, x, y, z));
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
        std::stringstream token(words[0]);
        int element_number;
        token >> element_number;
        std::stringstream t_type(words[1]);
        int tempvalue;
        t_type >> tempvalue;
        // assert(siconos::mechanics::fem::is_valid_element(tempvalue) &&
        //        "Unknown element type."); // --> check done in FiniteElementModel->init()
        FiniteElementType element_type = static_cast<FiniteElementType>(tempvalue);

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

        std::vector<std::shared_ptr<MeshVertex>> vertices_e;
        for (decltype(words.size()) k = 3 + number_of_tags; k < words.size(); k++) {
          decltype(vertices.size()) node_number;
          std::stringstream token(words[k]);
          token >> node_number;
          //        std::cout << "node_number" << node_number << std::endl;
          size_t v;
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
            std::make_shared<MeshElement>(element_number, element_type, vertices_e, tags));
        // elements.back()->display();
      }
    }
  }

  in.close();

  return std::make_shared<Mesh>(m, vertices, elements, physical_entities);
}

void siconos::mechanics::fem::writeMeshforPython(const Mesh& mesh,
                                                 std::string output_file_name) {
  std::filesystem::create_directories("outputs");
  std::filesystem::path filepath =
      std::filesystem::path("outputs") / ("mesh_" + output_file_name + ".py");

  std::ofstream outfile(filepath);
  if (!outfile.is_open())
    throw std::runtime_error("Cannot open file mesh_" + output_file_name + ".py");

  outfile << "coord = []\n";
  for (const auto& v : mesh.vertices()) {
    outfile << "coord.append([" << v->x() << "," << v->y() << "])\n";
  }

  outfile << "triangle = []\n";
  for (auto e : mesh.elements()) {
    outfile << "triangle.append([";
    for (auto v : e->vertices()) {
      outfile << v->num() << ",";
    }
    outfile << "])\n";
  }
  outfile.close();
}

std::string siconos::mechanics::fem::prepareWriteForPython(
    std::string basename, std::string fieldname, const std::vector<std::string>& fields) {
  std::string filename = basename + "_" + fieldname + ".py";
  std::filesystem::create_directories("outputs");

  std::filesystem::path filepath = std::filesystem::path("outputs") / filename;

  std::cout << "Output " << fieldname << " for python post-processing in ./outputs/"
            << filename << std::endl;
  std::ofstream outfile(filepath);
  if (!outfile.is_open()) throw std::runtime_error("Cannot open file " + filename);
  outfile << "import numpy as np\n";
  for (auto& str : fields) outfile << str + "=[]\n";
  outfile.close();
  return filename;
}

std::string siconos::mechanics::fem::prepareWriteDisplacementforPython(std::string basename) {
  return prepareWriteForPython(basename, "displacement", {"x", "y", "z"});
}

std::string siconos::mechanics::fem::prepareWriteTensorforPython(std::string basename,
                                                                 std::string tensorName) {
  return prepareWriteForPython(basename, "tensor",
                               {tensorName + "_xx", tensorName + "_yy", tensorName + "_xy"});
}

std::string siconos::mechanics::fem::prepareWriteBeamTensorforPython(std::string basename,
                                                                     std::string tensorName) {
  return prepareWriteForPython(
      basename, "tensor",
      {tensorName + "_tension", tensorName + "_bending1", tensorName + "_bending2"});
}

void siconos::mechanics::fem::writeDisplacementforPython(
    const Mesh& mesh, const FiniteElementModel& femodel,
    const siconos::algebra::SiconosVector& x, std::string filename) {
  std::filesystem::create_directories("outputs");

  std::filesystem::path filepath = std::filesystem::path("outputs") / filename;
  std::ofstream outfile(filepath, std::ios::app);
  if (!outfile.is_open()) throw std::runtime_error("Cannot open file " + filename);
  outfile.precision(15);
  outfile.setf(std::ios::scientific);

  std::vector<std::string> fields = {"x", "y", "z"};
  size_t pos = 0;
  for (auto& str : fields) {
    outfile << str + ".append(np.array([";
    if (pos == 2 && mesh.dim() < 3) {
      for (auto v : mesh.vertices()) {
        outfile << 0.0 << ", ";
      }
    } else {
      for (auto v : mesh.vertices()) {
        auto node = femodel.vertexToNode(v);
        double value = 0.0;
        if (node) {
          auto idx = node->global_dof_index()[pos];
          value = x(idx);
        }
        outfile << value << ", ";
      }
    }
    outfile << "]))\n\n";
    pos += 1;
  }
  outfile.close();
}

void siconos::mechanics::fem::writeTensorforPython(const FiniteElementModel& femodel,
                                                   const siconos::algebra::SiconosVector& x,
                                                   std::string filename,
                                                   std::string tensorName) {
  std::filesystem::path filepath = std::filesystem::path("outputs") / filename;
  std::ofstream outfile(filepath, std::ios::app);
  if (!outfile.is_open()) throw std::runtime_error("Cannot open file " + filename);
  outfile.precision(15);
  outfile.setf(std::ios::scientific);

  std::vector<std::string> fields = {tensorName + "_xx", tensorName + "_yy",
                                     tensorName + "_xy"};
  int pos = 0;
  for (auto& str : fields) {
    siconos::algebra::Index elem_cnt = 0;
    outfile << str + ".append(np.array([";
    for (size_t elem_cnt = 0; elem_cnt < femodel.elements().size(); elem_cnt++) {
      auto count = siconos::algebra::to_index(3 * elem_cnt);
      outfile << x(count + pos) << ", ";
    }
    std::cout << " elem_cnt:" << elem_cnt << std::endl;
    outfile << "]))\n\n";
    pos += 1;
  }
  outfile.close();
}

void siconos::mechanics::fem::prepareWriteBeamPositionforSOFA(std::string filename) {
  std::cout << "Output displacement for SOFA post-processing in " << filename << std::endl;

  FILE* foutput = fopen(filename.c_str(), "w");
  fclose(foutput);
}

void siconos::mechanics::fem::prepareWriteBlockPositionforSOFA(std::string filename) {
  std::cout << "Output displacement for SOFA post-processing in " << filename << std::endl;

  FILE* foutput = fopen(filename.c_str(), "w");
  fclose(foutput);
}

void siconos::mechanics::fem::writeBeamPositionforSOFA(
    const Mesh& mesh, const FiniteElementModel& femodel,
    const siconos::algebra::SiconosVector& positions, std::string filename, double t) {
  FILE* foutput = fopen(filename.c_str(), "a");
  fprintf(foutput, "T= %f\n", t);
  fprintf(foutput, "X=");
  for (const auto& v : mesh.vertices()) {
    std::shared_ptr<FENode> n = femodel.vertexToNode(v);
    double value = 0.0, valueqx, valueqy, valueqz;
    if (n) {
      auto idx = n->global_dof_index()[0];
      auto idy = n->global_dof_index()[1];

      value = positions(idx);
      fprintf(foutput, " %e", value + v->x());

      value = positions(idy);
      fprintf(foutput, " %e", value + v->y());

      if (n->global_dof_index().size() == 3) {
        // 2D case
        fprintf(foutput, " %e", 0.0);

        auto idtheta = n->global_dof_index()[2];
        // idqx= (*n->global_dof_index())[3];
        valueqx = positions(idtheta);
        // idqy= (*n->global_dof_index())[4];
        // valueqy =positions(idqy);
        // idqz= (*n->global_dof_index())[5];
        // valueqz =positions(idqz);

        double quat[4];

        double c3 = cos((valueqx + v->z()) / 2);
        double c1 = cos(0);
        double c2 = cos(0);

        double s3 = sin((valueqx + v->z()) / 2);
        double s1 = sin(0);
        double s2 = sin(0);

        quat[0] = s1 * c2 * c3 - c1 * s2 * s3;
        quat[1] = c1 * s2 * c3 + s1 * c2 * s3;
        quat[2] = c1 * c2 * s3 - s1 * s2 * c3;
        quat[3] = c1 * c2 * c3 + s1 * s2 * s3;

        fprintf(foutput, " %e", quat[0]);
        fprintf(foutput, " %e", quat[1]);
        fprintf(foutput, " %e", quat[2]);
        fprintf(foutput, " %e", quat[3]);
      } else {
        auto idz = n->global_dof_index()[2];
        value = positions(idz);
        fprintf(foutput, " %e", value + v->z());

        auto idqx = n->global_dof_index()[3];
        valueqx = positions(idqx);
        auto idqy = n->global_dof_index()[4];
        valueqy = positions(idqy);
        auto idqz = n->global_dof_index()[5];
        valueqz = positions(idqz);

        siconos::algebra::SiconosVector3 angle;
        angle << valueqx, valueqy, valueqz;
        auto quat = siconos::geometry::quaternionFromRotationVector(angle);

        fprintf(foutput, " %e", quat(4));
        fprintf(foutput, " %e", quat(5));
        fprintf(foutput, " %e", quat(6));
        fprintf(foutput, " %e", quat(3));
      }
    }
  }
  fprintf(foutput, "\n");
  fclose(foutput);
}

void siconos::mechanics::fem::writeBlockPositionforSOFA(
    const siconos::algebra::SiconosVector& positions, std::string filename, double t) {
  FILE* foutput = fopen(filename.c_str(), "a");
  fprintf(foutput, "T= %f\n", t);
  fprintf(foutput, "X=");
  fprintf(foutput, " %e", positions(0));
  fprintf(foutput, " %e", positions(1));
  // 2D case
  fprintf(foutput, " %e", 0.0);

  fprintf(foutput, "\n");
  fclose(foutput);
}
