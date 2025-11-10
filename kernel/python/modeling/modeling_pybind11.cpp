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
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <Eigen/Sparse>
#include <cstddef>

#include "DynamicalSystem.hpp"
#include "Interaction.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "NonSmoothLaw.hpp"
#include "Relation.hpp"
#include "SiconosMatrix.hpp"
// #include "eigen2python_pybind11.hpp"

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

void wrap_dynamical_systems(py::module_& m);
void wrap_nonsmoothlaws(py::module_& m);
void wrap_relations(py::module_& m);

std::vector<std::shared_ptr<siconos::modeling::Interaction>> interactions(
    std::shared_ptr<siconos::graphs::InteractionsGraph> dsg) {
  std::vector<std::shared_ptr<siconos::modeling::Interaction>> r =
      std::vector<std::shared_ptr<siconos::modeling::Interaction>>();
  siconos::graphs::InteractionsGraph::VIterator vi, viend;
  for (boost::tie(vi, viend) = dsg->vertices(); vi != viend; ++vi) {
    r.push_back(dsg->bundle(*vi));
  };
  return r;
};

siconos::algebra::SiconosSparseMatrix scale_sparse(siconos::algebra::SiconosSparseMatrix& m,
                                                   double scalar) {
  m *= scalar;
  return m;
}

py::object eigen_to_csc(siconos::algebra::SiconosSparseMatrix& M) {
  M.makeCompressed();
  auto nrows = M.rows();
  auto ncols = M.cols();
  auto nnz = M.nonZeros();

  // Création des arrays numpy sans copie
  auto data = py::array_t<double>(nnz, M.valuePtr());
  auto indices = py::array_t<siconos::algebra::SparseIndex>(nnz, M.innerIndexPtr());
  auto indptr = py::array_t<siconos::algebra::SparseIndex>(ncols + 1, M.outerIndexPtr());

  py::module scipy_sparse = py::module::import("scipy.sparse");
  py::object csc = scipy_sparse.attr("csc_array");

  return scipy_sparse.attr("csc_array")(py::make_tuple(data, indices, indptr),
                                        py::make_tuple(nrows, ncols));
}
siconos::algebra::SiconosSparseMatrix csc_to_eigen(py::object csc) {
  auto data =
      csc.attr("data").cast<py::array_t<double, py::array::c_style | py::array::forcecast>>();
  auto indices = csc.attr("indices")
                     .cast<py::array_t<siconos::algebra::SparseIndex,
                                       py::array::c_style | py::array::forcecast>>();
  auto indptr = csc.attr("indptr")
                    .cast<py::array_t<siconos::algebra::SparseIndex,
                                      py::array::c_style | py::array::forcecast>>();
  auto shape = csc.attr("shape").cast<std::pair<int64_t, int64_t>>();

  siconos::algebra::SparseIndex rows = shape.first, cols = shape.second;
  siconos::algebra::SparseIndex nnz = static_cast<siconos::algebra::SparseIndex>(data.size());

  auto data_ptr = data.mutable_data();
  auto idx_ptr = indices.mutable_data();
  auto indptr_ptr = indptr.mutable_data();

  std::vector<Eigen::Triplet<double, siconos::algebra::SparseIndex>> triplets;
  triplets.reserve(nnz);

  for (siconos::algebra::SparseIndex j = 0; j < cols; ++j) {
    auto start = indptr_ptr[j];
    auto end = indptr_ptr[j + 1];
    for (auto k = start; k < end; ++k) {
      triplets.emplace_back(idx_ptr[k], j, data_ptr[k]);
    }
  }
  siconos::algebra::SiconosSparseMatrix M(rows, cols);
  M.setFromTriplets(triplets.begin(), triplets.end());

  return M;  // RVO
}

PYBIND11_MODULE(modeling, m) {
  // py::module_ simulation_module = py::module_::import("siconos.simulation");

  // Optional docstring
  m.doc() = "Siconos modeling library";

  // m.def("scale_sparse", &scale_sparse, "Scale a sparse matrix by a scalar");
  // py::class_<siconos::modeling::LDS>(m, "LDS")
  //     .def(py::init<int>())
  //     .def("compute", &siconos::modeling::LDS::compute)
  //     .def("mass", &siconos::modeling::LDS::mass,
  //     py::return_value_policy::reference_internal) .def("reset",
  //     &siconos::modeling::LDS::reset) .def("check", &siconos::modeling::LDS::check) .def(
  //         "setCompute2",
  //         [](siconos::modeling::LDS& self, py::function f) {
  //           // Catch Python function and create a complient std::function
  //           self.setCompute([f](siconos::algebra::SiconosSparseMatrix& result) { f(result);
  //           });
  //         },
  //         "How to compute mass matrix")

  //     //  .def("setCompute",
  //     //       [](siconos::modeling::LDS& self, py::function pyfunc) {
  //     //         self._func = [pyfunc](siconos::algebra::SiconosSparseMatrix& M) {
  //     //           py::gil_scoped_acquire gil;
  //     //           py::object pyM = eigen_to_csc(M);
  //     //           pyM = pyfunc(pyM);
  //     //           M = siconos::pybind11_utils::csc_to_eigen(pyM);
  //     //         };
  //     //       })
  //     .def(
  //         "setComputeRuntime",
  //         [](siconos::modeling::LDS& self, py::function pyfunc) {
  //           self._func = [pyfunc](siconos::algebra::SiconosSparseMatrix& M) {
  //             pybind11::gil_scoped_acquire gil;
  //             pyfunc(M);
  //           };
  //           // self.setCompute(func);
  //         },
  //         py::keep_alive<1, 2>())
  //     .def("setComputeRuntime2",
  //          [](siconos::modeling::LDS& self, py::function pyfunc) {
  //            auto func = [pyfunc](siconos::algebra::SiconosSparseMatrix M)
  //                -> siconos::algebra::SiconosSparseMatrix {
  //              pybind11::gil_scoped_acquire gil;
  //              py::object result = pyfunc(std::move(M));
  //              return result.cast<siconos::algebra::SiconosSparseMatrix>();
  //            };
  //            self.setCompute2(func);
  //          })

  //     .def("setComputeRuntime4",
  //          [](siconos::modeling::LDS& self, py::function pyfunc) {
  //            auto func = [pyfunc](siconos::algebra::SiconosSparseMatrix& M)
  //                -> siconos::algebra::SiconosSparseMatrix {
  //              pybind11::gil_scoped_acquire gil;
  //              py::object result = pyfunc(M);
  //              return result.cast<siconos::algebra::SiconosSparseMatrix>();
  //            };
  //            self.setCompute4(func);
  //          })
  //     .def("setComputeRuntime3", [](siconos::modeling::LDS& self, py::function pyfunc) {
  //       auto func = [pyfunc](siconos::algebra::SiconosSparseMatrix& M)
  //           -> siconos::algebra::SiconosSparseMatrix {
  //         pybind11::gil_scoped_acquire gil;
  //         py::object result = pyfunc(M);
  //         return result.cast<siconos::algebra::SiconosSparseMatrix>();
  //       };
  //       self.setCompute3(func);
  //     });

  wrap_dynamical_systems(m);
  wrap_nonsmoothlaws(m);
  wrap_relations(m);

  // CLASSES with no Derived classes
  py::class_<siconos::modeling::Interaction, std::shared_ptr<siconos::modeling::Interaction>>(
      m, "Interaction")
      .def(py::init<std::shared_ptr<siconos::modeling::NonSmoothLaw>,
                    std::shared_ptr<siconos::modeling::Relation>>())
      .def("lambda_python", &siconos::modeling::Interaction::lambda_python,
           py::return_value_policy::reference_internal)
      .def("y", &siconos::modeling::Interaction::y_python,
           py::return_value_policy::reference_internal)
      .def("relation", &siconos::modeling::Interaction::relation,
           py::return_value_policy::reference_internal);

  m.def("interactions", &interactions, py::arg("graph"),
        "Return a list of Interaction objects from an InteractionsGraph");

  py::class_<siconos::modeling::NonSmoothDynamicalSystem,
             std::shared_ptr<siconos::modeling::NonSmoothDynamicalSystem>>(
      m, "NonSmoothDynamicalSystem")
      .def(py::init<double, double>())
      .def("insertDynamicalSystem",
           &siconos::modeling::NonSmoothDynamicalSystem::insertDynamicalSystem)
      .def("getNumberOfDS",
           &siconos::modeling::NonSmoothDynamicalSystem::getNumberOfDS,
	   "get the number of DS contained in the NonSmoothDynamicalSystem")
      .def("dynamicalSystem",
           &siconos::modeling::NonSmoothDynamicalSystem::dynamicalSystem,
	   "get the Ds number nb")
      .def("removeDynamicalSystem",
           &siconos::modeling::NonSmoothDynamicalSystem::removeDynamicalSystem,
	   py::arg("ds"),
	   "remove the Ds")

      .def("dynamicalSystemsVector",
           &siconos::modeling::NonSmoothDynamicalSystem::dynamicalSystemsVector,
	   "get the array of all DSs contained in the NonSmoothDynamicalSystem")

      .def("link", &siconos::modeling::NonSmoothDynamicalSystem::link,
           "link an interaction to two dynamical systems", py::arg("inter"), py::arg("ds1"),
           py::arg("ds2") = std::shared_ptr<siconos::modeling::DynamicalSystem>())

      .def("setTitle", &siconos::modeling::NonSmoothDynamicalSystem::setTitle, "set DS title")

      .def("topology", &siconos::modeling::NonSmoothDynamicalSystem::topology,
           "display the topology of the system")

      .def("setName",
           py::overload_cast<std::shared_ptr<siconos::modeling::DynamicalSystem>,
                             const std::string&>(
               &siconos::modeling::NonSmoothDynamicalSystem::setName),
           "set DS name")
      .def("setName",
           py::overload_cast<std::shared_ptr<siconos::modeling::Interaction>,
                             const std::string&>(
               &siconos::modeling::NonSmoothDynamicalSystem::setName),
           "set Interaction name")
      .def("displayDynamicalSystems",
           &siconos::modeling::NonSmoothDynamicalSystem::displayDynamicalSystems,
           "Print all dynamical systems infos")
      .def("__repr__", [](const siconos::modeling::NonSmoothDynamicalSystem& a) {
        a.display();
        return "\n";
      });
}
