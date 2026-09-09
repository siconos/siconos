/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2026 INRIA.
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

#include <pybind11/eigen.h>
#include <pybind11/native_enum.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include "FemTools.hpp"
#include "FiniteElementLinearTIDS.hpp"
#include "FiniteElementModel.hpp"
#include "LagrangianSparseDS.hpp"
#include "LagrangianSparseLinearTIDS.hpp"
#include "Material.hpp"
#include "Mesh.hpp"
#include "MeshUtils.hpp"

namespace py = pybind11;

PYBIND11_MODULE(_fem, m) {
  py::module_ modeling_module = py::module_::import("siconos.modeling");

  m.doc() = "Siconos mechanics.fem module";

  m.def("createMeshFromGMSH2", &siconos::mechanics::fem::createMeshFromGMSH2, py::arg("input"),
        py::arg("is_filename") = true);

  py::class_<siconos::mechanics::fem::Material, py::smart_holder> material(m, "Material");
  material
      .def(py::init<double, double, double, double>(),  // parameterized constructor
           py::arg("massDensity"), py::arg("E"), py::arg("nu"), py::arg("radius") = 1.0)
      .def("massDensity", &siconos::mechanics::fem::Material::massDensity)
      .def("elasticYoungModulus", &siconos::mechanics::fem::Material::elasticYoungModulus)
      .def("Poisson_s_ratio", &siconos::mechanics::fem::Material::Poisson_s_ratio)
      .def("analysisType2D", &siconos::mechanics::fem::Material::analysisType2D)
      .def("thickness", &siconos::mechanics::fem::Material::thickness)
      .def("shearModulus", &siconos::mechanics::fem::Material::shearModulus)
      .def("momentOfInertia", &siconos::mechanics::fem::Material::momentOfInertia)
      .def("crossSectionArea", &siconos::mechanics::fem::Material::crossSectionArea)
      .def("secondMomentOfArea", &siconos::mechanics::fem::Material::secondMomentOfArea);

  material.attr("Steel") = siconos::mechanics::fem::Steel

      ;

  py::native_enum<siconos::mechanics::fem::AnalysisType2D>(material, "AnalysisType2D",
                                                           "enum.Enum")
      // should be upper case, not PEP 8 compliant
      .value("plane_strain", siconos::mechanics::fem::AnalysisType2D::plane_strain)
      .value("plane_stress", siconos::mechanics::fem::AnalysisType2D::plane_stress)
      .value("axisymmetric", siconos::mechanics::fem::AnalysisType2D::axysymmetric)
      .export_values()
      .finalize();

  py::class_<siconos::mechanics::fem::FiniteElementLinearTIDS,
             siconos::modeling::LagrangianSparseLinearTIDS, py::smart_holder>(
      m, "FiniteElementLinearTIDS")
      .def(py::init<std::shared_ptr<siconos::mechanics::fem::Mesh>,
                    const std::map<int, const siconos::mechanics::fem::Material>&>(),
           py::arg("mesh"), py::arg("materials"))
      .def("FEModel", &siconos::mechanics::fem::FiniteElementLinearTIDS::FEModel)
      .def("applyDirichletBoundaryConditions",
           &siconos::mechanics::fem::FiniteElementLinearTIDS::applyDirichletBoundaryConditions,
           py::arg("physical_entity_tag"), py::arg("node_dof_index"),
           py::arg("imposedVelocity") = 0.0)
      .def("applyNodalForces",
           &siconos::mechanics::fem::FiniteElementLinearTIDS::applyNodalForces,
           py::arg("physical_entity_tag"), py::arg("nodal_forces"))
      .def("elasticPotentialEnergy",
           &siconos::mechanics::fem::FiniteElementLinearTIDS::elasticPotentialEnergy)
      .def("mass", &siconos::modeling::LagrangianSparseDS::mass_py)
      .def("stiffnessMatrix",
           &siconos::modeling::LagrangianSparseLinearTIDS::stiffnessMatrix_py)
      .def("fext_vector",
           [](siconos::mechanics::fem::FiniteElementLinearTIDS& self) {
             // Return a copy of the vector
             return siconos::algebra::SiconosVector(self.fext());
           })
      .def("q",
           [](siconos::mechanics::fem::FiniteElementLinearTIDS& self) { return *(self.q()); })
      .def("velocity",
           [](siconos::mechanics::fem::FiniteElementLinearTIDS& self) {
             return *(self.velocity());
           })
      .def("display", &siconos::mechanics::fem::FiniteElementLinearTIDS::display,
           py::arg("brief") = true);

  py::enum_<siconos::mechanics::fem::MeshTags>(m, "MeshTags")
      .value("bulk_material", siconos::mechanics::fem::MeshTags::bulk_material)
      .value("boundary_conditions", siconos::mechanics::fem::MeshTags::boundary_conditions)
      .value("applied_forces", siconos::mechanics::fem::MeshTags::applied_forces)
      .export_values();

  m.attr("default_tags") = siconos::mechanics::fem::default_tags;

  py::class_<siconos::mechanics::fem::FiniteElementModel, py::smart_holder>(
      m, "FiniteElementModel")
      .def("mesh", &siconos::mechanics::fem::FiniteElementModel::mesh)  // method
      .def("vertices",
           [](const siconos::mechanics::fem::FiniteElementModel& self) {
             // Convert span to vector for pybind11
             auto span = self.mesh()->vertices();
             return std::vector<std::shared_ptr<siconos::mechanics::fem::MeshVertex>>(
                 span.begin(), span.end());
           })
      .def("elements",
           [](const siconos::mechanics::fem::FiniteElementModel& self) {
             auto span = self.mesh()->elements();
             return std::vector<std::shared_ptr<siconos::mechanics::fem::MeshElement>>(
                 span.begin(), span.end());
           })
      .def("vertexToNode", &siconos::mechanics::fem::FiniteElementModel::vertexToNode)
      .def("nodes",
           [](const siconos::mechanics::fem::FiniteElementModel& self) {
             auto span = self.nodes();
             return std::vector<std::shared_ptr<siconos::mechanics::fem::FENode>>(span.begin(),
                                                                                  span.end());
           })
      .def("display", &siconos::mechanics::fem::FiniteElementModel::display,
           py::arg("brief") = true)
      .def("getContactSegments",
           &siconos::mechanics::fem::FiniteElementModel::getContactSegments,
           py::arg("contact_tag"))
      .def("getContactGlobalIndices",
           [](const siconos::mechanics::fem::FiniteElementModel& self,
              py::array_t<double> coords) {
             auto buf = coords.request();
             if (buf.ndim != 2 || buf.shape[1] != 2)
               throw std::runtime_error("contact_coords must be (N,2)");
             Eigen::Map<const Eigen::MatrixXd> mat(static_cast<double*>(buf.ptr), buf.shape[0],
                                                   2);
             return self.getContactGlobalIndices(mat);
           });

  py::class_<siconos::mechanics::fem::Mesh, py::smart_holder>(m, "Mesh")
      .def("dim", &siconos::mechanics::fem::Mesh::dim)
      .def("vertices",
           [](const siconos::mechanics::fem::Mesh& self) {
             auto span = self.vertices();
             return std::vector<std::shared_ptr<siconos::mechanics::fem::MeshVertex>>(
                 span.begin(), span.end());
           })
      .def("elements",
           [](const siconos::mechanics::fem::Mesh& self) {
             auto span = self.elements();
             return std::vector<std::shared_ptr<siconos::mechanics::fem::MeshElement>>(
                 span.begin(), span.end());
           })
      .def("physical_entities", &siconos::mechanics::fem::Mesh::physical_entities)
      .def("display", &siconos::mechanics::fem::Mesh::display, py::arg("brief") = true);

  py::class_<siconos::mechanics::fem::MeshVertex, py::smart_holder>(m, "MeshVertex")
      .def("num", &siconos::mechanics::fem::MeshVertex::num)
      .def("x", &siconos::mechanics::fem::MeshVertex::x)
      .def("y", &siconos::mechanics::fem::MeshVertex::y)
      .def("z", &siconos::mechanics::fem::MeshVertex::z)
      .def("display", &siconos::mechanics::fem::MeshVertex::display);

  py::class_<siconos::mechanics::fem::FENode, py::smart_holder>(m, "FENode")
      .def("num", &siconos::mechanics::fem::FENode::num)
      .def("global_dof_index",
           [](const siconos::mechanics::fem::FENode& self) {
             return std::vector<siconos::algebra::Index>(self.global_dof_index().begin(),
                                                         self.global_dof_index().end());
           })
      .def("x", &siconos::mechanics::fem::FENode::x)
      .def("y", &siconos::mechanics::fem::FENode::y)
      .def("z", &siconos::mechanics::fem::FENode::z)
      .def("display", &siconos::mechanics::fem::FENode::display);
}
