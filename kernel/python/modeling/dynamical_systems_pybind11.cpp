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

#include <Eigen/Sparse>
#include <cstddef>
#include <sstream>
// #include <functional>
// #include <memory>
// #include <span>

#include "FirstOrderLinearDS.hpp"
#include "FunctionTypes.hpp"
#include "LagrangianDS.hpp"
#include "LagrangianLinearDiagonalDS.hpp"
#include "LagrangianLinearTIDS.hpp"
#include "LagrangianSparseDS.hpp"
#include "LagrangianSparseLinearTIDS.hpp"
#include "NewtonEulerDS.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"
#include "eigen2python_pybind11.hpp"

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

void wrap_dynamical_systems(py::module_& m) {
  m.doc() = "Siconos modeling library";

  py::class_<siconos::modeling::DynamicalSystem,
             std::shared_ptr<siconos::modeling::DynamicalSystem>>(m, "DynamicalSystem")
      .def("x", &siconos::modeling::DynamicalSystem::x_python,
           py::return_value_policy::reference_internal)
      .def("r", &siconos::modeling::DynamicalSystem::r_python,
           py::return_value_policy::reference_internal)
      .def_property_readonly("dimension", &siconos::modeling::DynamicalSystem::dimension,
                             "dimension of the system")

      .def("resetCount", &siconos::modeling::DynamicalSystem::resetCount)
      .def("resetToInitialState", &siconos::modeling::DynamicalSystem::resetToInitialState)
      .def("swapInMemory", &siconos::modeling::DynamicalSystem::swapInMemory)
      .def("setNumber", &siconos::modeling::DynamicalSystem::setNumber)
      .def("number", &siconos::modeling::DynamicalSystem::number)
      .def("display", &siconos::modeling::DynamicalSystem::display)
      .def("__str__",
           [](const siconos::modeling::DynamicalSystem& self) {
             std::ostringstream buffer;
             py::scoped_ostream_redirect redirect(std::cout,
                                                  py::module_::import("sys").attr("stdout"));
             self.display(true);
             return buffer.str();
           })
      .def("__repr__", [](const siconos::modeling::DynamicalSystem& self) {
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
           py::arg("x0"))
      .def("setConstantMMatrixAlias",
           &siconos::modeling::FirstOrderNonLinearDS::setConstantMMatrix,
           py::keep_alive<1, 2>(), "To define a constant mass matrix");

  py::class_<siconos::modeling::FirstOrderLinearDS,
             std::shared_ptr<siconos::modeling::FirstOrderLinearDS>,
             siconos::modeling::FirstOrderNonLinearDS>(m, "FirstOrderLinearDS")
      .def(py::init<Eigen::Ref<siconos::algebra::SiconosVector>&>(), py::keep_alive<1, 2>(),
           py::arg("x0"))
      .def(py::init<Eigen::Ref<siconos::algebra::SiconosVector>&,
                    Eigen::Ref<siconos::algebra::SiconosDenseMatrix>&,
                    Eigen::Ref<siconos::algebra::SiconosVector>&>(),
           py::keep_alive<1, 2>(), py::keep_alive<1, 3>(), py::keep_alive<1, 4>(),
           py::arg("x0"), py::arg("A"), py::arg("b"))
      .def("setConstantA", &siconos::modeling::FirstOrderLinearDS::setConstantA,
           py::keep_alive<1, 2>(), "To define a constant A operator")
      .def("setConstantbVector", &siconos::modeling::FirstOrderLinearDS::setConstantbVector,
           py::keep_alive<1, 2>(), "To define a constant b operator")
      .def(
          "setComputebVectorFunction",
          [](siconos::modeling::FirstOrderLinearDS& self, py::function f) {
            // Catch Python function and create a complient std::function
            self.setComputebVectorFunction(
                [f](double val, Eigen::Ref<siconos::algebra::MapVectorType> result) {
                  f(val, result);  // Call python func with a memory view ...
                });
          },
          "How to compute b(t)");

  // ============================== SECOND ORDER DS ==============================

  py::class_<siconos::modeling::SecondOrderDS,
             std::shared_ptr<siconos::modeling::SecondOrderDS>,
             siconos::modeling::DynamicalSystem>(m, "SecondOrderDS")
      .def("p", &siconos::modeling::SecondOrderDS::p_python,
           py::return_value_policy::reference_internal);

  py::class_<siconos::modeling::LagrangianSparseDS,
             std::shared_ptr<siconos::modeling::LagrangianSparseDS>,
             siconos::modeling::SecondOrderDS>(m, "LagrangianSparseDS", py::dynamic_attr())
      .def(py::init<Eigen::Ref<siconos::algebra::SiconosVector>,
                    Eigen::Ref<siconos::algebra::SiconosVector>>(),
           py::keep_alive<1, 2>(),  // keep python object (np array arguments) memory alive
                                    // as long as object is referenced
           py::keep_alive<1, 3>(), py::arg("q0"), py::arg("v0"))
      .def("q", &siconos::modeling::LagrangianSparseDS::q_python,
           py::return_value_policy::reference_internal)
      .def("velocity", &siconos::modeling::LagrangianSparseDS::velocity_python,
           py::return_value_policy::reference_internal)
      .def("setConstantFext", &siconos::modeling::LagrangianSparseDS::setConstantFext,
           py::keep_alive<1, 2>(), "To define a constant external forces vector")
      .def(
          "setComputeFextFunction",
          [](siconos::modeling::LagrangianSparseDS& self, py::function f) {
            // Catch Python function and create a complient std::function
            self.setComputeFextFunction(
                [f](double val, Eigen::Ref<siconos::algebra::MapVectorType> result) {
                  f(val, result);  // Call python func with a memory view ...
                });
          },
          "How to compute external forces")

      .def("computeFext", &siconos::modeling::LagrangianSparseDS::computeFext,
           "compute external forces")

      .def("fext", &siconos::modeling::LagrangianSparseDS::fext,
           "current values of external forces")

      .def("setConstantMassCopy",
           [](siconos::modeling::LagrangianSparseDS& self, py::object csc) {
             py::object py_self = py::cast(self);
             if (py::hasattr(py_self, "_mass_ref")) py_self.attr("_mass_ref") = py::none();
             siconos::algebra::SiconosSparseMatrix M =
                 siconos::pybind11_utils::csc_to_eigen(csc);
             self.setConstantMassCopy(std::move(M));  // move => one allocation, no second copy
           })
      .def("setConstantMassAlias",
           [](siconos::modeling::LagrangianSparseDS& self, py::object csc) {
             auto map_ptr = siconos::pybind11_utils::csc_to_eigen_map(csc);
             self.setConstantMassAlias(*map_ptr);
             py::object py_self = py::cast(self, py::return_value_policy::reference);
             py_self.attr("_mass_ref") = csc;
           })
      .def_property(
          "mass_view",  // Read only
          [](const siconos::modeling::LagrangianSparseDS& self) -> py::object {
            if (!self.hasMass()) {
              return py::none();
            }
            py::object py_self = py::cast(self);
            // First case: memory has been set using an alias with a python/scipy matrix
            if (py::hasattr(py_self, "_mass_ref") && !py_self.attr("_mass_ref").is_none()) {
              // We choose to avoid this case :
              // - Mass set with an alias to a scipy mat.
              // - Read only mode can't be achieve
              // Don't want to mislead the user.
              throw std::runtime_error(
                  "Mass matrix has been set with an alias to a given scipy matrix. Please use "
                  "mass_alias function or the original scipy matrix to access to the mass "
                  "attribute.");
            }
            // Second case: memory comes from C++ side
            else {
              const auto& M = self.mass_py();
              // return siconos::pybind11_utils::eigensparse_to_scipy(M, py::cast(self));
              return siconos::pybind11_utils::make_readonly_csc_array(M, py::cast(self));
            }
          },
          nullptr, "Read-only view on the mass (scipy.sparse.csc_array).")

      .def_property(
          "mass_alias",
          [](const siconos::modeling::LagrangianSparseDS& self) -> py::object {
            if (!self.hasMass()) {
              return py::none();
            }

            py::object py_self = py::cast(self);
            // First case: memory has been set using an alias with a python/scipy matrix
            if (py::hasattr(py_self, "_mass_ref") && !py_self.attr("_mass_ref").is_none()) {
              py::object sp = py::module_::import("scipy.sparse");
              auto csc = py_self.attr("_mass_ref");
              if (!(py::isinstance(csc, sp.attr("csc_array")))) {
                throw std::runtime_error(
                    "_mass_ref is corrupted: expected scipy.sparse.csc_array");
              }

              return csc;  // In that case we get directly the original scipy matrix
              // Warning: the result won't be read-only!!
            }
            // Second case: memory comes from C++ side
            else {
              throw std::runtime_error(
                  "Mass matrix has not been set with an alias to a given scipy matrix. Please "
                  "use "
                  "mass_view function to access (read-only) to the mass attribute.");
            }
          },
          nullptr, "View (shared memory) on the mass (scipy.sparse.csc_array).")

      .def(
          "setComputeMassFunction",
          [](siconos::modeling::LagrangianSparseDS& self, py::function pyfunc) {
            auto func = [pyfunc](const Eigen::Ref<const siconos::algebra::SiconosVector>& q,
                                 siconos::algebra::SiconosSparseMatrix& M)
                -> siconos::algebra::SiconosSparseMatrix {
              pybind11::gil_scoped_acquire gil;
              py::object result = pyfunc(q, M);
              return result.cast<siconos::algebra::SiconosSparseMatrix>();
            };
            self.setComputeMassFunction(func);
          },
          "How to compute mass matrix")

      .def("computeMass", &siconos::modeling::LagrangianSparseDS::computeMass)
      .def("hasMass", &siconos::modeling::LagrangianSparseDS::hasMass)
      .def("hasConstantMass", &siconos::modeling::LagrangianSparseDS::hasConstantMass)
      .def("setConstantJacobianFintOver_qCopy",
           [](siconos::modeling::LagrangianSparseDS& self, py::object csc) {
             py::object py_self = py::cast(self);
             if (py::hasattr(py_self, "_jacobianFintOver_q_ref"))
               py_self.attr("_jacobianFintOver_q_ref") = py::none();

             siconos::algebra::SiconosSparseMatrix M =
                 siconos::pybind11_utils::csc_to_eigen(csc);
             self.setConstantJacobianFintOver_qCopy(
                 std::move(M));  // move => one allocation, no second copy
           })
      .def("setConstantJacobianFintOver_qAlias",
           [](siconos::modeling::LagrangianSparseDS& self, py::object csc) {
             auto map_ptr = siconos::pybind11_utils::csc_to_eigen_map(csc);
             self.setConstantJacobianFintOver_qAlias(*map_ptr);
             py::object py_self = py::cast(self, py::return_value_policy::reference);
             py_self.attr("_jacobianFintOver_q_ref") = csc;
           })

      .def_property(
          "jacobianFintOver_q_view",  // Read only
          [](const siconos::modeling::LagrangianSparseDS& self) -> py::object {
            if (!self.hasJacobianFintOver_q()) {
              return py::none();
            }
            py::object py_self = py::cast(self);
            // First case: memory has been set using an alias with a python/scipy matrix
            if (py::hasattr(py_self, "_jacobianFintOver_q_ref") &&
                !py_self.attr("_jacobianFintOver_q_ref").is_none()) {
              // We choose to avoid this case :
              // - JacobianFintOver_q set with an alias to a scipy mat.
              // - Read only mode can't be achieve
              // Don't want to mislead the user.
              throw std::runtime_error(
                  "JacobianFintOver_q matrix has been set with an alias to a given scipy "
                  "matrix. Please use "
                  "jacobianFintOver_q_alias function or the original scipy matrix to access "
                  "to the jacobianFintOver_q "
                  "attribute.");
            }
            // Second case: memory comes from C++ side
            else {  //(self.owned_jacobianFintOver_q_mat_) {
              const auto& M = self.jacobianFintOver_q_py();
              return siconos::pybind11_utils::make_readonly_csc_array(M, py_self);
            }
          },
          nullptr, "Read-only view on the jacobianFintOver_q (scipy.sparse.csc_array).")

      .def_property(
          "jacobianFintOver_q_alias",
          [](const siconos::modeling::LagrangianSparseDS& self) -> py::object {
            if (!self.hasJacobianFintOver_q()) {
              return py::none();
            }

            py::object py_self = py::cast(self);
            // First case: memory has been set using an alias with a python/scipy matrix
            if (py::hasattr(py_self, "_jacobianFintOver_q_ref") &&
                !py_self.attr("_jacobianFintOver_q_ref").is_none()) {
              py::object sp = py::module_::import("scipy.sparse");
              auto csc = py_self.attr("_jacobianFintOver_q_ref");
              if (!(py::isinstance(csc, sp.attr("csc_array")))) {
                throw std::runtime_error(
                    "_jacobianFintOver_q_ref is corrupted: expected scipy.sparse.csc_array");
              }

              return csc;  // In that case we get directly the original scipy matrix
              // Warning: the result won't be read-only!!
            }
            // Second case: memory comes from C++ side
            else {
              throw std::runtime_error(
                  "JacobianFintOver_q matrix has not been set with an alias to a given scipy "
                  "matrix. Please "
                  "use "
                  "jacobianFintOver_q_view function to access (read-only) to the "
                  "jacobianFintOver_q attribute.");
            }
          },
          nullptr, "View (shared memory) on the jacobianFintOver_q (scipy.sparse.csc_array).")

      .def("setComputeJacobianFintOver_qFunction",
           [](siconos::modeling::LagrangianSparseDS& self, py::function pyfunc) {
             auto func = [pyfunc](
                             const Eigen::Ref<const siconos::algebra::SiconosVector>& velocity,
                             const Eigen::Ref<const siconos::algebra::SiconosVector>& position,
                             double current_time, siconos::algebra::SiconosSparseMatrix& M)
                 -> siconos::algebra::SiconosSparseMatrix {
               pybind11::gil_scoped_acquire gil;
               py::object result = pyfunc(velocity, position, current_time, M);
               return result.cast<siconos::algebra::SiconosSparseMatrix>();
             };
             self.setComputeJacobianFintOver_qFunction(func);
           })

      .def("computeJacobianFintOver_q",
           &siconos::modeling::LagrangianSparseDS::computeJacobianFintOver_q)
      .def("hasJacobianFintOver_q",
           &siconos::modeling::LagrangianSparseDS::hasJacobianFintOver_q)
      .def("hasConstantJacobianTotalForcesOver_q",
           &siconos::modeling::LagrangianSparseDS::hasConstantJacobianTotalForcesOver_q)
      .def("hasConstantJacobianTotalForcesOver_velocity",
           &siconos::modeling::LagrangianSparseDS::hasConstantJacobianTotalForcesOver_velocity)

      .def("setConstantJacobianFintOver_velocityCopy",
           [](siconos::modeling::LagrangianSparseDS& self, py::object csc) {
             py::object py_self = py::cast(self);
             if (py::hasattr(py_self, "_jacobianFintOver_velocity_ref"))
               py_self.attr("_jacobianFintOver_velocity_ref") = py::none();

             siconos::algebra::SiconosSparseMatrix M =
                 siconos::pybind11_utils::csc_to_eigen(csc);
             self.setConstantJacobianFintOver_velocityCopy(
                 std::move(M));  // move => one allocation, no second copy
           })
      .def("setConstantJacobianFintOver_velocityAlias",
           [](siconos::modeling::LagrangianSparseDS& self, py::object csc) {
             auto map_ptr = siconos::pybind11_utils::csc_to_eigen_map(csc);
             self.setConstantJacobianFintOver_velocityAlias(*map_ptr);
             py::object py_self = py::cast(self, py::return_value_policy::reference);
             py_self.attr("_jacobianFintOver_velocity_ref") = csc;
           })

      .def_property(
          "jacobianFintOver_velocity_view",  // Read only
          [](const siconos::modeling::LagrangianSparseDS& self) -> py::object {
            if (!self.hasJacobianFintOver_velocity()) {
              return py::none();
            }
            py::object py_self = py::cast(self);
            // First case: memory has been set using an alias with a python/scipy matrix
            if (py::hasattr(py_self, "_jacobianFintOver_velocity_ref") &&
                !py_self.attr("_jacobianFintOver_velocity_ref").is_none()) {
              // We choose to avoid this case :
              // - JacobianFintOver_velocity set with an alias to a scipy mat.
              // - Read only mode can't be achieve
              // Don't want to mislead the user.
              throw std::runtime_error(
                  "JacobianFintOver_velocity matrix has been set with an alias to a given "
                  "scipy "
                  "matrix. Please use "
                  "jacobianFintOver_velocity_alias function or the original scipy matrix to "
                  "access "
                  "to the jacobianFintOver_velocity "
                  "attribute.");
            }
            // Second case: memory comes from C++ side
            else {  //(self.owned_jacobianFintOver_velocity_mat_) {
              const auto& M = self.jacobianFintOver_velocity_py();
              return siconos::pybind11_utils::make_readonly_csc_array(M, py_self);
            }
          },
          nullptr, "Read-only view on the jacobianFintOver_velocity (scipy.sparse.csc_array).")

      .def_property(
          "jacobianFintOver_velocity_alias",
          [](const siconos::modeling::LagrangianSparseDS& self) -> py::object {
            if (!self.hasJacobianFintOver_velocity()) {
              return py::none();
            }

            py::object py_self = py::cast(self);
            // First case: memory has been set using an alias with a python/scipy matrix
            if (py::hasattr(py_self, "_jacobianFintOver_velocity_ref") &&
                !py_self.attr("_jacobianFintOver_velocity_ref").is_none()) {
              py::object sp = py::module_::import("scipy.sparse");
              auto csc = py_self.attr("_jacobianFintOver_velocity_ref");
              if (!(py::isinstance(csc, sp.attr("csc_array")))) {
                throw std::runtime_error(
                    "_jacobianFintOver_velocity_ref is corrupted: expected "
                    "scipy.sparse.csc_array");
              }

              return csc;  // In that case we get directly the original scipy matrix
              // Warning: the result won't be read-only!!
            }
            // Second case: memory comes from C++ side
            else {
              throw std::runtime_error(
                  "JacobianFintOver_velocity matrix has not been set with an alias to a given "
                  "scipy "
                  "matrix. Please "
                  "use "
                  "jacobianFintOver_velocity_view function to access (read-only) to the "
                  "jacobianFintOver_velocity attribute.");
            }
          },
          nullptr,
          "View (shared memory) on the jacobianFintOver_velocity (scipy.sparse.csc_array).")

      .def("setComputeJacobianFintOver_velocityFunction",
           [](siconos::modeling::LagrangianSparseDS& self, py::function pyfunc) {
             auto func = [pyfunc](
                             const Eigen::Ref<const siconos::algebra::SiconosVector>& velocity,
                             const Eigen::Ref<const siconos::algebra::SiconosVector>& position,
                             double current_time, siconos::algebra::SiconosSparseMatrix& M)
                 -> siconos::algebra::SiconosSparseMatrix {
               pybind11::gil_scoped_acquire gil;
               py::object result = pyfunc(velocity, position, current_time, M);
               return result.cast<siconos::algebra::SiconosSparseMatrix>();
             };
             self.setComputeJacobianFintOver_velocityFunction(func);
           })

      .def("computeJacobianFintOver_velocity",
           &siconos::modeling::LagrangianSparseDS::computeJacobianFintOver_velocity)
      .def("hasJacobianFintOver_velocity",
           &siconos::modeling::LagrangianSparseDS::hasJacobianFintOver_velocity)

      .def("setConstantJacobianFgyrOver_qCopy",
           [](siconos::modeling::LagrangianSparseDS& self, py::object csc) {
             py::object py_self = py::cast(self);
             if (py::hasattr(py_self, "_jacobianFgyrOver_q_ref"))
               py_self.attr("_jacobianFgyrOver_q_ref") = py::none();

             siconos::algebra::SiconosSparseMatrix M =
                 siconos::pybind11_utils::csc_to_eigen(csc);
             self.setConstantJacobianFgyrOver_qCopy(
                 std::move(M));  // move => one allocation, no second copy
           })
      .def("setConstantJacobianFgyrOver_qAlias",
           [](siconos::modeling::LagrangianSparseDS& self, py::object csc) {
             auto map_ptr = siconos::pybind11_utils::csc_to_eigen_map(csc);
             self.setConstantJacobianFgyrOver_qAlias(*map_ptr);
             py::object py_self = py::cast(self, py::return_value_policy::reference);
             py_self.attr("_jacobianFgyrOver_q_ref") = csc;
           })

      .def_property(
          "jacobianFgyrOver_q_view",  // Read only
          [](const siconos::modeling::LagrangianSparseDS& self) -> py::object {
            if (!self.hasJacobianFgyrOver_q()) {
              return py::none();
            }
            py::object py_self = py::cast(self);
            // First case: memory has been set using an alias with a python/scipy matrix
            if (py::hasattr(py_self, "_jacobianFgyrOver_q_ref") &&
                !py_self.attr("_jacobianFgyrOver_q_ref").is_none()) {
              // We choose to avoid this case :
              // - JacobianFgyrOver_q set with an alias to a scipy mat.
              // - Read only mode can't be achieve
              // Don't want to mislead the user.
              throw std::runtime_error(
                  "JacobianFgyrOver_q matrix has been set with an alias to a given scipy "
                  "matrix. Please use "
                  "jacobianFgyrOver_q_alias function or the original scipy matrix to access "
                  "to the jacobianFgyrOver_q "
                  "attribute.");
            }
            // Second case: memory comes from C++ side
            else {  //(self.owned_jacobianFgyrOver_q_mat_) {
              const auto& M = self.jacobianFgyrOver_q_py();
              return siconos::pybind11_utils::make_readonly_csc_array(M, py_self);
            }
          },
          nullptr, "Read-only view on the jacobianFgyrOver_q (scipy.sparse.csc_array).")

      .def_property(
          "jacobianFgyrOver_q_alias",
          [](const siconos::modeling::LagrangianSparseDS& self) -> py::object {
            if (!self.hasJacobianFgyrOver_q()) {
              return py::none();
            }

            py::object py_self = py::cast(self);
            // First case: memory has been set using an alias with a python/scipy matrix
            if (py::hasattr(py_self, "_jacobianFgyrOver_q_ref") &&
                !py_self.attr("_jacobianFgyrOver_q_ref").is_none()) {
              py::object sp = py::module_::import("scipy.sparse");
              auto csc = py_self.attr("_jacobianFgyrOver_q_ref");
              if (!(py::isinstance(csc, sp.attr("csc_array")))) {
                throw std::runtime_error(
                    "_jacobianFgyrOver_q_ref is corrupted: expected scipy.sparse.csc_array");
              }

              return csc;  // In that case we get directly the original scipy matrix
              // Warning: the result won't be read-only!!
            }
            // Second case: memory comes from C++ side
            else {
              throw std::runtime_error(
                  "JacobianFgyrOver_q matrix has not been set with an alias to a given scipy "
                  "matrix. Please "
                  "use "
                  "jacobianFgyrOver_q_view function to access (read-only) to the "
                  "jacobianFgyrOver_q attribute.");
            }
          },
          nullptr, "View (shared memory) on the jacobianFgyrOver_q (scipy.sparse.csc_array).")

      .def("setComputeJacobianFgyrOver_qFunction",
           [](siconos::modeling::LagrangianSparseDS& self, py::function pyfunc) {
             auto func = [pyfunc](
                             const Eigen::Ref<const siconos::algebra::SiconosVector>& velocity,
                             const Eigen::Ref<const siconos::algebra::SiconosVector>& position,
                             siconos::algebra::SiconosSparseMatrix& M)
                 -> siconos::algebra::SiconosSparseMatrix {
               pybind11::gil_scoped_acquire gil;
               py::object result = pyfunc(velocity, position, M);
               return result.cast<siconos::algebra::SiconosSparseMatrix>();
             };
             self.setComputeJacobianFgyrOver_qFunction(func);
           })

      .def("computeJacobianFgyrOver_q",
           &siconos::modeling::LagrangianSparseDS::computeJacobianFgyrOver_q)
      .def("hasJacobianFgyrOver_q",
           &siconos::modeling::LagrangianSparseDS::hasJacobianFgyrOver_q)

      .def("setConstantJacobianFgyrOver_velocityCopy",
           [](siconos::modeling::LagrangianSparseDS& self, py::object csc) {
             py::object py_self = py::cast(self);
             if (py::hasattr(py_self, "_jacobianFgyrOver_velocity_ref"))
               py_self.attr("_jacobianFgyrOver_velocity_ref") = py::none();

             siconos::algebra::SiconosSparseMatrix M =
                 siconos::pybind11_utils::csc_to_eigen(csc);
             self.setConstantJacobianFgyrOver_velocityCopy(
                 std::move(M));  // move => one allocation, no second copy
           })
      .def("setConstantJacobianFgyrOver_velocityAlias",
           [](siconos::modeling::LagrangianSparseDS& self, py::object csc) {
             auto map_ptr = siconos::pybind11_utils::csc_to_eigen_map(csc);
             self.setConstantJacobianFgyrOver_velocityAlias(*map_ptr);
             py::object py_self = py::cast(self, py::return_value_policy::reference);
             py_self.attr("_jacobianFgyrOver_velocity_ref") = csc;
           })

      .def_property(
          "jacobianFgyrOver_velocity_view",  // Read only
          [](const siconos::modeling::LagrangianSparseDS& self) -> py::object {
            if (!self.hasJacobianFgyrOver_velocity()) {
              return py::none();
            }
            py::object py_self = py::cast(self);
            // First case: memory has been set using an alias with a python/scipy matrix
            if (py::hasattr(py_self, "_jacobianFgyrOver_velocity_ref") &&
                !py_self.attr("_jacobianFgyrOver_velocity_ref").is_none()) {
              // We choose to avoid this case :
              // - JacobianFgyrOver_velocity set with an alias to a scipy mat.
              // - Read only mode can't be achieve
              // Don't want to mislead the user.
              throw std::runtime_error(
                  "JacobianFgyrOver_velocity matrix has been set with an alias to a given "
                  "scipy "
                  "matrix. Please use "
                  "jacobianFgyrOver_velocity_alias function or the original scipy matrix to "
                  "access "
                  "to the jacobianFgyrOver_velocity "
                  "attribute.");
            }
            // Second case: memory comes from C++ side
            else {  //(self.owned_jacobianFgyrOver_velocity_mat_) {
              const auto& M = self.jacobianFgyrOver_velocity_py();
              return siconos::pybind11_utils::make_readonly_csc_array(M, py_self);
            }
          },
          nullptr, "Read-only view on the jacobianFgyrOver_velocity (scipy.sparse.csc_array).")

      .def_property(
          "jacobianFgyrOver_velocity_alias",
          [](const siconos::modeling::LagrangianSparseDS& self) -> py::object {
            if (!self.hasJacobianFgyrOver_velocity()) {
              return py::none();
            }

            py::object py_self = py::cast(self);
            // First case: memory has been set using an alias with a python/scipy matrix
            if (py::hasattr(py_self, "_jacobianFgyrOver_velocity_ref") &&
                !py_self.attr("_jacobianFgyrOver_velocity_ref").is_none()) {
              py::object sp = py::module_::import("scipy.sparse");
              auto csc = py_self.attr("_jacobianFgyrOver_velocity_ref");
              if (!(py::isinstance(csc, sp.attr("csc_array")))) {
                throw std::runtime_error(
                    "_jacobianFgyrOver_velocity_ref is corrupted: expected "
                    "scipy.sparse.csc_array");
              }

              return csc;  // In that case we get directly the original scipy matrix
              // Warning: the result won't be read-only!!
            }
            // Second case: memory comes from C++ side
            else {
              throw std::runtime_error(
                  "JacobianFgyrOver_velocity matrix has not been set with an alias to a given "
                  "scipy "
                  "matrix. Please "
                  "use "
                  "jacobianFgyrOver_velocity_view function to access (read-only) to the "
                  "jacobianFgyrOver_velocity attribute.");
            }
          },
          nullptr,
          "View (shared memory) on the jacobianFgyrOver_velocity (scipy.sparse.csc_array).")

      .def("setComputeJacobianFgyrOver_velocityFunction",
           [](siconos::modeling::LagrangianSparseDS& self, py::function pyfunc) {
             auto func = [pyfunc](
                             const Eigen::Ref<const siconos::algebra::SiconosVector>& velocity,
                             const Eigen::Ref<const siconos::algebra::SiconosVector>& position,
                             siconos::algebra::SiconosSparseMatrix& M)
                 -> siconos::algebra::SiconosSparseMatrix {
               pybind11::gil_scoped_acquire gil;
               py::object result = pyfunc(velocity, position, M);
               return result.cast<siconos::algebra::SiconosSparseMatrix>();
             };
             self.setComputeJacobianFgyrOver_velocityFunction(func);
           })

      .def("computeJacobianFgyrOver_velocity",
           &siconos::modeling::LagrangianSparseDS::computeJacobianFgyrOver_velocity)
      .def("hasJacobianFgyrOver_velocity",
           &siconos::modeling::LagrangianSparseDS::hasJacobianFgyrOver_velocity);

  py::class_<siconos::modeling::LagrangianSparseLinearTIDS,
             std::shared_ptr<siconos::modeling::LagrangianSparseLinearTIDS>,
             siconos::modeling::LagrangianSparseDS>(m, "LagrangianSparseLinearTIDS")
      .def(py::init<Eigen::Ref<siconos::algebra::SiconosVector>,
                    Eigen::Ref<siconos::algebra::SiconosVector>,
                    const siconos::algebra::SiconosSparseMatrix&>(),
           py::keep_alive<1, 2>(),  // keep python object (np array arguments) memory
                                    // alive as long as object is referenced
           py::keep_alive<1, 3>(), py::arg("q0"), py::arg("v0"), py::arg("mass"),
           R"pbdoc(
              constructor from initial state and mass matrix only.
              warning: M will be copied into mass attribute.
         )pbdoc")

      .def(py::init<Eigen::Ref<siconos::algebra::SiconosVector>,
                    Eigen::Ref<siconos::algebra::SiconosVector>>(),
           py::keep_alive<1, 2>(),  // keep python object (np array arguments) memory
                                    // alive as long as object is referenced
           py::keep_alive<1, 3>(), py::arg("q0"), py::arg("v0"),
           R"pbdoc(constructor from initial state only.)pbdoc")
      .def("setStiffnessMatrixCopy",
           [](siconos::modeling::LagrangianSparseLinearTIDS& self, py::object csc) {
             py::object py_self = py::cast(self);
             if (py::hasattr(py_self, "_stiffness_ref"))
               py_self.attr("_stiffness_ref") = py::none();

             siconos::algebra::SiconosSparseMatrix K =
                 siconos::pybind11_utils::csc_to_eigen(csc);
             self.setStiffnessMatrixCopy(
                 std::move(K));  // move => one allocation, no second copy
           })

      .def("setStiffnessMatrixAlias",
           [](siconos::modeling::LagrangianSparseLinearTIDS& self, py::object csc) {
             auto map_ptr = siconos::pybind11_utils::csc_to_eigen_map(csc);
             self.setStiffnessMatrixAlias(*map_ptr);
             py::object py_self = py::cast(self, py::return_value_policy::reference);
             py_self.attr("_stiffness_ref") = csc;
           })

      .def(
          "setDampingMatrixCopy",
          [](siconos::modeling::LagrangianSparseLinearTIDS& self, py::object csc) {
            py::object py_self = py::cast(self);
            if (py::hasattr(py_self, "_damping_ref"))
              py_self.attr("_damping_ref") = py::none();

            siconos::algebra::SiconosSparseMatrix C =
                siconos::pybind11_utils::csc_to_eigen(csc);
            self.setDampingMatrixCopy(std::move(C));  // move => one allocation, no second copy
          })

      .def("setDampingMatrixAlias",
           [](siconos::modeling::LagrangianSparseLinearTIDS& self, py::object csc) {
             auto map_ptr = siconos::pybind11_utils::csc_to_eigen_map(csc);
             self.setDampingMatrixAlias(*map_ptr);
             py::object py_self = py::cast(self, py::return_value_policy::reference);
             py_self.attr("_damping_ref") = csc;
           })

      .def_property(
          "stiffness_view",  // Read only
          [](const siconos::modeling::LagrangianSparseLinearTIDS& self) -> py::object {
            if (!self.hasStiffnessMatrix()) {
              return py::none();
            }
            py::object py_self = py::cast(self);
            // First case: memory has been set using an alias with a python/scipy matrix
            if (py::hasattr(py_self, "_stiffness_ref") &&
                !py_self.attr("_stiffness_ref").is_none()) {
              // We choose to avoid this case :
              // - Matrix set with an alias to a scipy mat.
              // - Read only mode can't be achieve
              // Don't want to mislead the user.
              throw std::runtime_error(
                  "Stiffness matrix has been set with an alias to a given scipy matrix. "
                  "Please use "
                  "stiffness_view function or the original scipy matrix to access to the "
                  "stiffness "
                  "attribute.");
            }
            // Second case: memory comes from C++ side
            else {
              const auto& M = self.stiffnessMatrix_py();
              return siconos::pybind11_utils::make_readonly_csc_array(M, py_self);
            }
          },
          nullptr, "Read-only view on the stiffness matrix (scipy.sparse.csc_array).")

      .def_property(
          "stiffness_alias",
          [](const siconos::modeling::LagrangianSparseLinearTIDS& self) -> py::object {
            if (!self.hasStiffnessMatrix()) {
              return py::none();
            }

            py::object py_self = py::cast(self);
            // First case: memory has been set using an alias with a python/scipy matrix
            if (py::hasattr(py_self, "_stiffness_ref") &&
                !py_self.attr("_stiffness_ref").is_none()) {
              py::object sp = py::module_::import("scipy.sparse");
              auto csc = py_self.attr("_stiffness_ref");
              if (!(py::isinstance(csc, sp.attr("csc_array")))) {
                throw std::runtime_error(
                    "_stiffness_ref is corrupted: expected scipy.sparse.csc_array");
              }

              return csc;  // In that case we get directly the original scipy matrix
              // Warning: the result won't be read-only!!
            }
            // Second case: memory comes from C++ side
            else {
              throw std::runtime_error(
                  "Stiffness matrix has not been set with an alias to a given scipy matrix. "
                  "Please "
                  "use "
                  "stiffness_view function to access (read-only) to the stiffness attribute.");
            }
          },
          nullptr, "View (shared memory) on the stiffness (scipy.sparse.csc_array).")

      .def_property(
          "damping_view",  // Read only
          [](const siconos::modeling::LagrangianSparseLinearTIDS& self) -> py::object {
            if (!self.hasDampingMatrix()) {
              return py::none();
            }
            py::object py_self = py::cast(self);
            // First case: memory has been set using an alias with a python/scipy matrix
            if (py::hasattr(py_self, "_damping_ref") &&
                !py_self.attr("_damping_ref").is_none()) {
              // We choose to avoid this case :
              // - Matrix set with an alias to a scipy mat.
              // - Read only mode can't be achieve
              // Don't want to mislead the user.
              throw std::runtime_error(
                  "Damping matrix has been set with an alias to a given scipy matrix. "
                  "Please use "
                  "damping_view function or the original scipy matrix to access to the "
                  "damping "
                  "attribute.");
            }
            // Second case: memory comes from C++ side
            else {
              const auto& M = self.dampingMatrix_py();
              return siconos::pybind11_utils::make_readonly_csc_array(M, py_self);
            }
          },
          nullptr, "Read-only view on the damping matrix (scipy.sparse.csc_array).")

      .def_property(
          "damping_alias",
          [](const siconos::modeling::LagrangianSparseLinearTIDS& self) -> py::object {
            if (!self.hasDampingMatrix()) {
              return py::none();
            }

            py::object py_self = py::cast(self);
            // First case: memory has been set using an alias with a python/scipy matrix
            if (py::hasattr(py_self, "_damping_ref") &&
                !py_self.attr("_damping_ref").is_none()) {
              py::object sp = py::module_::import("scipy.sparse");
              auto csc = py_self.attr("_damping_ref");
              if (!(py::isinstance(csc, sp.attr("csc_array")))) {
                throw std::runtime_error(
                    "_damping_ref is corrupted: expected scipy.sparse.csc_array");
              }

              return csc;  // In that case we get directly the original scipy matrix
              // Warning: the result won't be read-only!!
            }
            // Second case: memory comes from C++ side
            else {
              throw std::runtime_error(
                  "Damping matrix has not been set with an alias to a given scipy matrix. "
                  "Please "
                  "use "
                  "damping_view function to access (read-only) to the damping attribute.");
            }
          },
          nullptr, "View (shared memory) on the damping (scipy.sparse.csc_array).");

  py::class_<siconos::modeling::LagrangianDS, std::shared_ptr<siconos::modeling::LagrangianDS>,
             siconos::modeling::SecondOrderDS>(m, "LagrangianDS")

      .def(py::init<Eigen::Ref<siconos::algebra::SiconosVector>,
                    Eigen::Ref<siconos::algebra::SiconosVector>>(),
           py::keep_alive<1, 2>(),  // keep python object (np array arguments) memory
                                    // alive as long as object is referenced
           py::keep_alive<1, 3>(), py::arg("q0"), py::arg("v0"))

      .def("q", &siconos::modeling::LagrangianDS::q_python,
           py::return_value_policy::reference_internal)
      .def("q0", &siconos::modeling::LagrangianDS::q0_python,
           py::return_value_policy::reference_internal)

      .def("velocity", &siconos::modeling::LagrangianDS::velocity_python,
           py::return_value_policy::reference_internal)
      .def("velocity0", &siconos::modeling::LagrangianDS::velocity0_python,
           py::return_value_policy::reference_internal)

      .def("setConstantFext", &siconos::modeling::LagrangianDS::setConstantFext,
           py::keep_alive<1, 2>(), "To define a constant external forces vector")

      .def(
          "setComputeFextFunction",
          [](siconos::modeling::LagrangianDS& self, py::function f) {
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
          [](siconos::modeling::LagrangianDS& self, py::function f) {
            // Catch Python function and create a complient std::function
            self.setComputeFintFunction(
                [f](const Eigen::Ref<const siconos::algebra::SiconosVector>& velocity,
                    const Eigen::Ref<const siconos::algebra::SiconosVector>& position,
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
          [](siconos::modeling::LagrangianDS& self, py::function f) {
            // Catch Python function and create a complient std::function
            self.setComputeFgyrFunction(
                [f](const Eigen::Ref<const siconos::algebra::SiconosVector>& velocity,
                    const Eigen::Ref<const siconos::algebra::SiconosVector>& position,
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

      .def("setConstantMassAlias", &siconos::modeling::LagrangianDS::setConstantMassAlias,
           py::keep_alive<1, 2>(), "To define a constant mass matrix")

      .def(
          "setComputeMassFunction",
          [](siconos::modeling::LagrangianDS& self, py::function f) {
            // Catch Python function and create a complient std::function
            self.setComputeMassFunction(
                [f](const Eigen::Ref<const siconos::algebra::SiconosVector>& q,
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
                    Eigen::Ref<siconos::algebra::SiconosDenseMatrix>>(),
           py::keep_alive<1, 2>(),  // keep python object (np array arguments) memory
                                    // alive as long as object is referenced
           py::keep_alive<1, 3>(), py::keep_alive<1, 4>(), py::arg("q0"), py::arg("v0"),
           py::arg("M"))

      .def("setStiffnessMatrix", &siconos::modeling::LagrangianLinearTIDS::setStiffnessMatrix,
           py::keep_alive<1, 2>(), "To define the stiffness matrix (constant)")
      .def("setDampingMatrix", &siconos::modeling::LagrangianLinearTIDS::setDampingMatrix,
           py::keep_alive<1, 2>(), "To define the damping matrix (constant)")
      .def_property_readonly("stiffnessMatrix",
                             &siconos::modeling::LagrangianLinearTIDS::stiffnessMatrix,
                             "stiffness matrix")
      .def_property_readonly("dampingMatrix",
                             &siconos::modeling::LagrangianLinearTIDS::dampingMatrix,
                             "damping matrix");

  py::class_<siconos::modeling::LagrangianLinearDiagonalDS,
             std::shared_ptr<siconos::modeling::LagrangianLinearDiagonalDS>,
             siconos::modeling::LagrangianDS>(m, "LagrangianLinearDiagonalDS")
      .def(py::init<Eigen::Ref<siconos::algebra::SiconosVector>,
                    Eigen::Ref<siconos::algebra::SiconosVector>,
                    Eigen::Ref<siconos::algebra::SiconosVector>>(),
           py::keep_alive<1, 2>(),  // keep python object (np array arguments) memory
                                    // alive as long as object is referenced
           py::keep_alive<1, 3>(), py::keep_alive<1, 4>(), py::arg("q0"), py::arg("v0"),
           py::arg("stiffness_diag"))

      .def_property_readonly("massMatrix",
                             &siconos::modeling::LagrangianLinearDiagonalDS::massMatrix,
                             "mass matrix")
      .def_property_readonly("stiffnessMatrix",
                             &siconos::modeling::LagrangianLinearDiagonalDS::stiffnessMatrix,
                             "stiffness matrix")
      .def_property_readonly("dampingMatrix",
                             &siconos::modeling::LagrangianLinearDiagonalDS::dampingMatrix,
                             "damping matrix");

  py::class_<siconos::modeling::NewtonEulerDS,
             std::shared_ptr<siconos::modeling::NewtonEulerDS>,
             siconos::modeling::SecondOrderDS>(m, "NewtonEulerDS")

      .def(py::init<Eigen::Ref<siconos::algebra::SiconosVector7>,
                    Eigen::Ref<siconos::algebra::SiconosVector6>, double,
                    Eigen::Ref<siconos::algebra::SiconosDenseMatrix>>(),
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
          [](siconos::modeling::NewtonEulerDS& self, py::function f) {
            // Catch Python function and create a complient std::function
            self.setComputeFextFunction(
                [f](double val, Eigen::Ref<siconos::algebra::MapVector3Type> result) {
                  f(val, result);  // Call python func with a memory view ...
                });
          },
          "How to compute external forces")
      .def("setConstantMext", &siconos::modeling::NewtonEulerDS::setConstantMext,
           py::keep_alive<1, 2>(), "To define a constant external torques vector")
      .def(
          "setComputeMextFunction",
          [](siconos::modeling::NewtonEulerDS& self, py::function f) {
            // Catch Python function and create a complient std::function
            self.setComputeMextFunction(
                [f](double val, Eigen::Ref<siconos::algebra::MapVector3Type> result) {
                  f(val, result);  // Call python func with a memory view ...
                });
          },
          "How to compute external torques")
      .def("setIsMextExpressedInInertialFrame",
           &siconos::modeling::NewtonEulerDS::setIsMextExpressedInInertialFrame)
      .def("setComputeJacobianMintOver_q_byFD",
           &siconos::modeling::NewtonEulerDS::setComputeJacobianMintOver_q_byFD)
      .def_property("scalarMass", &siconos::modeling::NewtonEulerDS::scalarMass,
                    &siconos::modeling::NewtonEulerDS::setScalarMass)
      .def("angularVelocity", &siconos::modeling::NewtonEulerDS::angularVelocity_view)
      .def("angularVelocityInBodyFrame",
           &siconos::modeling::NewtonEulerDS::angularVelocityInBodyFrame)
      .def_property_readonly("totalInertiaMatrix",
                             &siconos::modeling::NewtonEulerDS::totalInertiaMatrix,
                             "total Inertia matrix");
  ;
}
