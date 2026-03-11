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

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>

#include "PlasticityProblem.h"
#include "Plasticity_cst.h"
#include "NumericsMatrix.h"

namespace py = pybind11;

void wrap_plasticity(py::module_& m) {
  // Plasticity_DruckerPrager_model
  py::class_<Plasticity_DruckerPrager_model>(m, "Plasticity_DruckerPrager_model",
    "Drucker-Prager plasticity model parameters (eta, theta)")
    .def(py::init<>())
    .def_readwrite("eta", &Plasticity_DruckerPrager_model::eta, "Cone coefficients")
    .def_readwrite("theta", &Plasticity_DruckerPrager_model::theta, "Dilatancy coefficients");

  // PlasticityProblem - use a wrapper approach similar to friction contact
  py::class_<PlasticityProblem>(m, "PlasticityProblem",
    "Generic plasticity problem structure")
    .def(py::init<>())
    .def_readwrite("dimension", &PlasticityProblem::dimension, "Dimension of stress space")
    .def_readwrite("numberOfCones", &PlasticityProblem::numberOfCones, "Number of cones")
    .def_readwrite("M", &PlasticityProblem::M, "M matrix")
    .def_readwrite("q", &PlasticityProblem::q, "q vector")
    .def_readwrite("model", &PlasticityProblem::model, "Model parameters");

  // Factory function for creating PlasticityProblem with Drucker-Prager model
  m.def("PlasticityProblem_new", [](int dim, py::array_t<double> M, py::array_t<double> q,
                                     py::array_t<double> eta, py::array_t<double> theta) {
      auto* problem = new PlasticityProblem();
      problem->dimension = dim;
      problem->numberOfCones = eta.size();
      
      // Create M matrix from numpy array
      py::buffer_info M_info = M.request();
      problem->M = NM_create(NM_DENSE, M_info.shape[0], M_info.shape[1]);
      problem->M->matrix0 = (double*)malloc(M_info.shape[0] * M_info.shape[1] * sizeof(double));
      std::memcpy(problem->M->matrix0, M_info.ptr, M_info.shape[0] * M_info.shape[1] * sizeof(double));
      
      // Copy q vector
      py::buffer_info q_info = q.request();
      problem->q = (double*)malloc(q_info.shape[0] * sizeof(double));
      std::memcpy(problem->q, q_info.ptr, q_info.shape[0] * sizeof(double));
      
      // Create model and copy eta/theta
      problem->model = (Plasticity_DruckerPrager_model*)malloc(sizeof(Plasticity_DruckerPrager_model));
      py::buffer_info eta_info = eta.request();
      py::buffer_info theta_info = theta.request();
      problem->model->eta = (double*)malloc(eta_info.shape[0] * sizeof(double));
      problem->model->theta = (double*)malloc(theta_info.shape[0] * sizeof(double));
      std::memcpy(problem->model->eta, eta_info.ptr, eta_info.shape[0] * sizeof(double));
      std::memcpy(problem->model->theta, theta_info.ptr, theta_info.shape[0] * sizeof(double));
      
      return problem;
    }, py::arg("dim"), py::arg("M"), py::arg("q"), py::arg("eta"), py::arg("theta"),
       "Create a PlasticityProblem with Drucker-Prager model",
       py::return_value_policy::take_ownership);

  // Backward compatibility alias
  m.attr("Plasticity2DProblem") = m.attr("PlasticityProblem");
  
  // Backward compatibility factory function - same as PlasticityProblem_new
  m.def("Plasticity2DProblem_new", [](int dim, py::array_t<double> M, py::array_t<double> q,
                                       py::array_t<double> eta, py::array_t<double> theta) {
      auto* problem = new PlasticityProblem();
      problem->dimension = dim;
      problem->numberOfCones = eta.size();
      
      // Create M matrix from numpy array
      py::buffer_info M_info = M.request();
      problem->M = NM_create(NM_DENSE, M_info.shape[0], M_info.shape[1]);
      problem->M->matrix0 = (double*)malloc(M_info.shape[0] * M_info.shape[1] * sizeof(double));
      std::memcpy(problem->M->matrix0, M_info.ptr, M_info.shape[0] * M_info.shape[1] * sizeof(double));
      
      // Copy q vector
      py::buffer_info q_info = q.request();
      problem->q = (double*)malloc(q_info.shape[0] * sizeof(double));
      std::memcpy(problem->q, q_info.ptr, q_info.shape[0] * sizeof(double));
      
      // Create model and copy eta/theta
      problem->model = (Plasticity_DruckerPrager_model*)malloc(sizeof(Plasticity_DruckerPrager_model));
      py::buffer_info eta_info = eta.request();
      py::buffer_info theta_info = theta.request();
      problem->model->eta = (double*)malloc(eta_info.shape[0] * sizeof(double));
      problem->model->theta = (double*)malloc(theta_info.shape[0] * sizeof(double));
      std::memcpy(problem->model->eta, eta_info.ptr, eta_info.shape[0] * sizeof(double));
      std::memcpy(problem->model->theta, theta_info.ptr, theta_info.shape[0] * sizeof(double));
      
      return problem;
    }, py::arg("dim"), py::arg("M"), py::arg("q"), py::arg("eta"), py::arg("theta"),
       "Create a PlasticityProblem (backward compatibility)",
       py::return_value_policy::take_ownership);

  // Driver function (forward declaration - actual implementation is in the driver)
  // We need to declare it here to make it available in Python
  // The driver function signature is:
  // int plasticity_2d_driver(PlasticityProblem* problem, double* stress, double* plastic_strain_rate, SolverOptions* options);
  
  // Display functions
  m.def("plasticity_display", &plasticity_display,
        "Display plasticity problem",
        py::arg("problem"));

  // Backward compatibility function names
  m.def("mohrCoulomb2D_display", [](PlasticityProblem* problem) {
      plasticity_display(problem);
    }, "Display plasticity problem (backward compatibility)",
    py::arg("problem"));
  
  // Note: plasticity_2d_driver is not exposed in Python bindings because 
  // it requires numpy array handling for stress and plastic_strain_rate.
  // Users should use the C API or implement a custom wrapper.

  // Solver IDs enum
  py::enum_<PLASTICITY_SOLVER>(m, "PLASTICITY_SOLVER", "Plasticity solver IDs")
      .value("PLASTICITY_2D_NSGS", PLASTICITY_2D_NSGS, "Non-smooth Gauss-Seidel")
      .value("PLASTICITY_2D_ONECONE_NSN", PLASTICITY_2D_ONECONE_NSN, "Non-smooth Newton")
      .value("PLASTICITY_2D_ONECONE_NSN_GP", PLASTICITY_2D_ONECONE_NSN_GP, "Non-smooth Newton GP")
      .value("PLASTICITY_2D_ONECONE_NSN_GP_HYBRID", PLASTICITY_2D_ONECONE_NSN_GP_HYBRID, "Hybrid NSN")
      .value("PLASTICITY_2D_ONECONE_ProjectionOnCone", PLASTICITY_2D_ONECONE_ProjectionOnCone, "Projection on cone")
      .value("PLASTICITY_2D_ONECONE_ProjectionOnConeWithLocalIteration", 
             PLASTICITY_2D_ONECONE_ProjectionOnConeWithLocalIteration, "Projection with local iteration")
      .export_values();

  // Also export as module attributes for convenience
  m.attr("PLASTICITY_2D_NSGS") = (int)PLASTICITY_2D_NSGS;
  m.attr("PLASTICITY_2D_ONECONE_NSN") = (int)PLASTICITY_2D_ONECONE_NSN;
  m.attr("PLASTICITY_2D_ONECONE_NSN_GP") = (int)PLASTICITY_2D_ONECONE_NSN_GP;
  m.attr("PLASTICITY_2D_ONECONE_NSN_GP_HYBRID") = (int)PLASTICITY_2D_ONECONE_NSN_GP_HYBRID;
  m.attr("PLASTICITY_2D_ONECONE_ProjectionOnCone") = (int)PLASTICITY_2D_ONECONE_ProjectionOnCone;
  m.attr("PLASTICITY_2D_ONECONE_ProjectionOnConeWithLocalIteration") = 
      (int)PLASTICITY_2D_ONECONE_ProjectionOnConeWithLocalIteration;

  // Backward compatibility solver IDs
  m.attr("MOHR_COULOMB_2D_NSGS") = (int)PLASTICITY_2D_NSGS;
  m.attr("SICONOS_MC2D_NSGS") = (int)PLASTICITY_2D_NSGS;
}
