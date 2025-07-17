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

#include <pybind11/eigen.h>
#include <pybind11/iostream.h>
#include <pybind11/pybind11.h>

#include <sstream>
// #include <functional>
// #include <memory>
// #include <span>

#include "FirstOrderLinearDS.hpp"
#include "LagrangianDS.hpp"
#include "LagrangianLinearTIDS.hpp"
#include "LagrangianSparseDS.hpp"
#include "NewtonEulerDS.hpp"

// #include <pybind11/stl.h>  // Pour permettre la conversion entre std::vector et les objets
// Python comme les listes

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

void wrap_dynamical_systems(py::module_ &m) {
  m.doc() = "Siconos modeling library";

  py::class_<siconos::modeling::DynamicalSystem,
             std::shared_ptr<siconos::modeling::DynamicalSystem>>(m, "DynamicalSystem")
      .def("x", &siconos::modeling::DynamicalSystem::x_python,
           py::return_value_policy::reference_internal)
      .def("r", &siconos::modeling::DynamicalSystem::r_python,
           py::return_value_policy::reference_internal)
      .def("resetCount", &siconos::modeling::DynamicalSystem::resetCount)
      .def("setNumber", &siconos::modeling::DynamicalSystem::setNumber)
      .def("display", &siconos::modeling::DynamicalSystem::display)
      .def("__str__",
           [](const siconos::modeling::DynamicalSystem &self) {
             std::ostringstream buffer;
             py::scoped_ostream_redirect redirect(std::cout,
                                                  py::module_::import("sys").attr("stdout"));
             self.display(true);
             return buffer.str();
           })
      .def("__repr__", [](const siconos::modeling::DynamicalSystem &self) {
        std::ostringstream buffer;
        py::scoped_ostream_redirect redirect(std::cout,
                                             py::module_::import("sys").attr("stdout"));
        self.display(true);
        return buffer.str();
      });
  // ============================== FIRST ORDER DS ==============================

  py::class_<siconos::modeling::FirstOrderNonLinearDS,
             std::shared_ptr<siconos::modeling::FirstOrderNonLinearDS>,
             siconos::modeling::DynamicalSystem>(m, "FirstOrderNonLinearDS")
      .def(py::init<Eigen::Ref<siconos::algebra::SiconosVector>>(),
           py::keep_alive<1, 2>(),  // keep python object (np array arguments) memory alive
                                    // as long as object is referenced
           py::arg("x0"));

  py::class_<siconos::modeling::FirstOrderLinearDS,
             std::shared_ptr<siconos::modeling::FirstOrderLinearDS>,
             siconos::modeling::FirstOrderNonLinearDS>(m, "FirstOrderLinearDS")
      .def(py::init<Eigen::Ref<siconos::algebra::SiconosVector> &>(), py::keep_alive<1, 2>(),
           py::arg("x0"))
      .def(py::init<Eigen::Ref<siconos::algebra::SiconosVector> &,
                    Eigen::Ref<siconos::algebra::SiconosMatrix> &,
                    Eigen::Ref<siconos::algebra::SiconosVector> &>(),
           py::keep_alive<1, 2>(), py::keep_alive<1, 3>(), py::keep_alive<1, 4>(),
           py::arg("x0"), py::arg("A"), py::arg("b"))
      .def("setConstantA", &siconos::modeling::FirstOrderLinearDS::setConstantA,
           py::keep_alive<1, 2>(), "To define a constant A operator");

  // ============================== SECOND ORDER DS ==============================

  py::class_<siconos::modeling::SecondOrderDS,
             std::shared_ptr<siconos::modeling::SecondOrderDS>,
             siconos::modeling::DynamicalSystem>(m, "SecondOrderDS")
      .def("p", &siconos::modeling::SecondOrderDS::p_python,
           py::return_value_policy::reference_internal);
  py::class_<siconos::modeling::LagrangianSparseDS,
             std::shared_ptr<siconos::modeling::LagrangianSparseDS>,
             siconos::modeling::SecondOrderDS>(m, "LagrangianSparseDS")

      .def(py::init<Eigen::Ref<siconos::algebra::SiconosVector>,
                    Eigen::Ref<siconos::algebra::SiconosVector>>(),
           py::keep_alive<1, 2>(),  // keep python object (np array arguments) memory alive
                                    // as long as object is referenced
           py::keep_alive<1, 3>(), py::arg("q0"), py::arg("v0"))

      .def("q", &siconos::modeling::LagrangianSparseDS::q_python,
           py::return_value_policy::reference_internal)
      .def("velocity", &siconos::modeling::LagrangianSparseDS::velocity_python,
           py::return_value_policy::reference_internal)
      //  .def(
      //      "setConstantMass",
      //      [](siconos::modeling::LagrangianSparseDS &self,
      //         const siconos::algebra::SiconosSparseMatrix &mat) {
      //        self.setConstantMass(Eigen::Ref<const
      //        siconos::algebra::SiconosSparseMatrix>(mat));
      //      },
      //      py::arg("mass_matrix"))

      //  // Expose version shared_ptr (zéro copie, pour usage expert)
      //  .def(
      //      "setConstantMass",
      //      [](siconos::modeling::LagrangianSparseDS &self,
      //         std::shared_ptr<siconos::algebra::SiconosSparseMatrix> mat) {
      //        self.setConstantMass(mat);
      //      },
      //      py::arg("mass_matrix_shared"))
      //  .def(
      //      "setConstantMass",
      //      [](siconos::modeling::LagrangianSparseDS &self,
      //         const siconos::algebra::SiconosSparseMatrix &mat) {
      //        self.setConstantMass(Eigen::Ref<const
      //        siconos::algebra::SiconosSparseMatrix>(mat));
      //      },
      //      py::arg("mass_matrix"))

      .def("setConstantMass", &siconos::modeling::LagrangianSparseDS::setConstantMass,
           py::keep_alive<1, 2>(), "To define a constant mass operator")
      //           [](siconos::modeling::LagrangianSparseDS &self,
      //              const siconos::algebra::SiconosSparseMatrix &mat) {
      //     self.setConstantMass(std::make_shared<siconos::algebra::SiconosSparseMatrix>(mat));
      //           },
      //           py::arg("mass_matrix"))

      .def(
          "mass",
          [](const siconos::modeling::LagrangianSparseDS &self)
              -> const siconos::algebra::SiconosSparseMatrix & { return self.mass_py(); },
          py::return_value_policy::reference_internal)

      ;

  py::class_<siconos::modeling::LagrangianDS, std::shared_ptr<siconos::modeling::LagrangianDS>,
             siconos::modeling::SecondOrderDS>(m, "LagrangianDS")

      .def(py::init<Eigen::Ref<siconos::algebra::SiconosVector>,
                    Eigen::Ref<siconos::algebra::SiconosVector>>(),
           py::keep_alive<1, 2>(),  // keep python object (np array arguments) memory
                                    // alive as long as object is referenced
           py::keep_alive<1, 3>(), py::arg("q0"), py::arg("v0"))

      .def("q", &siconos::modeling::LagrangianDS::q_python,
           py::return_value_policy::reference_internal)
      .def("velocity", &siconos::modeling::LagrangianDS::velocity_python,
           py::return_value_policy::reference_internal)

      .def("setConstantFext", &siconos::modeling::LagrangianDS::setConstantFext,
           py::keep_alive<1, 2>(), "To define a constant external forces vector")

      .def(
          "setComputeFextFunction",
          [](siconos::modeling::LagrangianDS &self, py::function f) {
            // Catch Python function and create a complient std::function
            self.setComputeFextFunction(
                [f](double val, Eigen::Ref<siconos::algebra::MapVectorType> result) {
                  f(val, result);  // Call python func with a memory view ...
                });
          },
          "How to compute external forces")

      .def("computeFext", &siconos::modeling::LagrangianDS::computeFext,
           "compute external forces")

      .def("fext", &siconos::modeling::LagrangianDS::fext, "current values of external forces")

      .def(
          "setComputeFintFunction",
          [](siconos::modeling::LagrangianDS &self, py::function f) {
            // Catch Python function and create a complient std::function
            self.setComputeFintFunction(
                [f](const Eigen::Ref<const siconos::algebra::SiconosVector> &velocity,
                    const Eigen::Ref<const siconos::algebra::SiconosVector> &position,
                    double time, Eigen::Ref<siconos::algebra::MapVectorType> result) {
                  f(velocity, position, time,
                    result);  // Call python func with a memory view ...
                });
          },

          "How to compute internal forces")
      .def("computeFint", &siconos::modeling::LagrangianDS::computeFint,
           "compute internal forces")
      .def_property_readonly("fint", &siconos::modeling::LagrangianDS::fint,
                             "internal forces current values")

      .def(
          "setComputeFgyrFunction",
          [](siconos::modeling::LagrangianDS &self, py::function f) {
            // Catch Python function and create a complient std::function
            self.setComputeFgyrFunction(
                [f](const Eigen::Ref<const siconos::algebra::SiconosVector> &velocity,
                    const Eigen::Ref<const siconos::algebra::SiconosVector> &position,
                    Eigen::Ref<siconos::algebra::MapVectorType> result) {
                  f(velocity, position,
                    result);  // Call python func with a memory view ...
                });
          },

          "How to compute gyroscopic forces")
      .def("computeFgyr", &siconos::modeling::LagrangianDS::computeFgyr,
           "compute gyroscopic forces")
      .def_property_readonly("fgyr", &siconos::modeling::LagrangianDS::fgyr,
                             "gyroscopic forces current values")

      .def("computeTotalForces", &siconos::modeling::LagrangianDS::computeTotalForces,
           "compute total forces")
      .def_property_readonly("totalForces", &siconos::modeling::LagrangianDS::totalForces,
                             "internal forces current values")
      .def("computeJacobianTotalForcesOver_q",
           &siconos::modeling::LagrangianDS::computeJacobianTotalForcesOver_q,
           "compute total forces")
      .def_property_readonly("jacobianTotalForcesOver_q",
                             &siconos::modeling::LagrangianDS::jacobianTotalForcesOver_q,
                             "jacobian of the internal forces over q")
      .def("computeJacobianTotalForcesOver_velocity",
           &siconos::modeling::LagrangianDS::computeJacobianTotalForcesOver_velocity,
           "compute total forces")
      .def_property_readonly(
          "jacobianTotalForcesOver_velocity",
          &siconos::modeling::LagrangianDS::jacobianTotalForcesOver_velocity,
          "jacobian of the internal forces over velocity")

      .def("setConstantMass", &siconos::modeling::LagrangianDS::setConstantMass,
           py::keep_alive<1, 2>(), "To define a constant mass matrix")

      .def(
          "setComputeMassFunction",
          [](siconos::modeling::LagrangianDS &self, py::function f) {
            // Catch Python function and create a complient std::function
            self.setComputeMassFunction(
                [f](const Eigen::Ref<const siconos::algebra::SiconosVector> &q,
                    Eigen::Ref<siconos::algebra::MapType> result) {
                  f(q, result);  // Call python func with a memory view ...
                });
          },
          "How to compute mass matrix")

      .def("computeMass", &siconos::modeling::LagrangianDS::computeMass,
           "compute total forces")

      .def_property_readonly("mass", &siconos::modeling::LagrangianDS::mass, "mass matrix");

  py::class_<siconos::modeling::LagrangianLinearTIDS,
             std::shared_ptr<siconos::modeling::LagrangianLinearTIDS>,
             siconos::modeling::LagrangianDS>(m, "LagrangianLinearTIDS")
      .def(py::init<Eigen::Ref<siconos::algebra::SiconosVector>,
                    Eigen::Ref<siconos::algebra::SiconosVector>,
                    Eigen::Ref<siconos::algebra::SiconosMatrix>>(),
           py::keep_alive<1, 2>(),  // keep python object (np array arguments) memory
                                    // alive as long as object is referenced
           py::keep_alive<1, 3>(), py::keep_alive<1, 4>(), py::arg("q0"), py::arg("v0"),
           py::arg("M"))

      .def("setStiffnessMatrix", &siconos::modeling::LagrangianLinearTIDS::setStiffnessMatrix,
           py::keep_alive<1, 2>(), "To define the stiffness matrix (constant)")
      .def("setDampingMatrix", &siconos::modeling::LagrangianLinearTIDS::setDampingMatrix,
           py::keep_alive<1, 2>(), "To define the damping matrix (constant)");

  py::class_<siconos::modeling::NewtonEulerDS,
             std::shared_ptr<siconos::modeling::NewtonEulerDS>,
             siconos::modeling::SecondOrderDS>(m, "NewtonEulerDS")

      .def(py::init<Eigen::Ref<siconos::algebra::SiconosVector>,
                    Eigen::Ref<siconos::algebra::SiconosVector>, double,
                    Eigen::Ref<siconos::algebra::SiconosMatrix>>(),
           py::keep_alive<1, 2>(),  // keep python object (np array arguments) memory
                                    // alive as long as object is referenced
           py::keep_alive<1, 3>(), py::keep_alive<1, 5>(), py::arg("q0"), py::arg("twist0"),
           py::arg("mass"), py::arg("inertia"))

      .def("q", &siconos::modeling::NewtonEulerDS::q_python,
           py::return_value_policy::reference_internal)
      //  .def("velocity", &siconos::modeling::NewtonEulerDS::twist_python,
      //       py::return_value_policy::reference_internal)
      .def("twist", &siconos::modeling::NewtonEulerDS::twist_python,
           py::return_value_policy::reference_internal)

      .def("setConstantFext", &siconos::modeling::NewtonEulerDS::setConstantFext,
           py::keep_alive<1, 2>(), "To define a constant external forces vector")
      .def(
          "setComputeFextFunction",
          [](siconos::modeling::NewtonEulerDS &self, py::function f) {
            // Catch Python function and create a complient std::function
            self.setComputeFextFunction(
                [f](double val, Eigen::Ref<siconos::algebra::MapVectorType> result) {
                  f(val, result);  // Call python func with a memory view ...
                });
          },
          "How to compute external forces")
      .def("setConstantMext", &siconos::modeling::NewtonEulerDS::setConstantMext,
           py::keep_alive<1, 2>(), "To define a constant external torques vector")
      .def(
          "setComputeMextFunction",
          [](siconos::modeling::NewtonEulerDS &self, py::function f) {
            // Catch Python function and create a complient std::function
            self.setComputeMextFunction(
                [f](double val, Eigen::Ref<siconos::algebra::MapVectorType> result) {
                  f(val, result);  // Call python func with a memory view ...
                });
          },
          "How to compute external torques")
      .def("setIsMextExpressedInInertialFrame",
           &siconos::modeling::NewtonEulerDS::setIsMextExpressedInInertialFrame)
      .def_property("scalarMass", &siconos::modeling::NewtonEulerDS::scalarMass,
                    &siconos::modeling::NewtonEulerDS::setScalarMass)
      .def("angularVelocity", &siconos::modeling::NewtonEulerDS::angularVelocity_view)
      .def("angularVelocityInBodyFrame",
           &siconos::modeling::NewtonEulerDS::angularVelocityInBodyFrame);
}
