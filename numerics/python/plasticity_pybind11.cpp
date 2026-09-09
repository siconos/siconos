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

#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>

#include "NumericsMatrix.h"
#include "PlasticityProblem.h"
#include "Plasticity_options.h"

namespace py = pybind11;

void wrap_plasticity(py::module_& m, py::module_& params, py::module_& solver_ids) {
  // PlasticityModelType enum
  py::enum_<PlasticityModelType>(m, "PlasticityModelType", "Type of plasticity model")
      .value("UNKNOWN", PLASTICITY_MODEL_UNKNOWN, "Unknown model")
      .value("DRUCKER_PRAGER", PLASTICITY_MODEL_DRUCKER_PRAGER, "Drucker-Prager model")
      .value("VON_MISES", PLASTICITY_MODEL_VON_MISES, "Von Mises model")
      .export_values();

  // Plasticity_DruckerPrager_model
  py::class_<Plasticity_DruckerPrager_model, py::smart_holder>(
      m, "Plasticity_DruckerPrager_model",
      "Drucker-Prager plasticity model parameters (eta, theta)")
      .def(py::init<>())
      .def_readwrite("eta", &Plasticity_DruckerPrager_model::eta, "Cone coefficients")
      .def_readwrite("theta", &Plasticity_DruckerPrager_model::theta,
                     "Dilatancy coefficients");

  // Plasticity_VonMises_model
  py::class_<Plasticity_VonMises_model, py::smart_holder>(
      m, "Plasticity_VonMises_model", "Von Mises plasticity model parameters (yield stress)")
      .def(py::init<>())
      .def_readwrite("sigma_y", &Plasticity_VonMises_model::sigma_y, "Yield stress");

  // PlasticityProblem - use a wrapper approach similar to friction contact
  py::class_<PlasticityProblem, py::smart_holder>(m, "PlasticityProblem",
                                                  "Generic plasticity problem structure")
      .def(py::init<>())
      .def_readwrite("dimension", &PlasticityProblem::dimension, "Dimension of stress space")
      .def_readwrite("numberOfCones", &PlasticityProblem::numberOfCones, "Number of cones")
      .def_readwrite("M", &PlasticityProblem::M, "M matrix")
      .def_readwrite("q", &PlasticityProblem::q, "q vector")
      .def_readwrite("model_type", &PlasticityProblem::model_type, "Model type enum")
      .def_readwrite("model", &PlasticityProblem::model, "Model parameters union");

  // Factory function for creating PlasticityProblem with Drucker-Prager model
  m.def(
      "PlasticityProblem_new",
      [](int dim, py::array_t<double> M, py::array_t<double> q, py::array_t<double> eta,
         py::array_t<double> theta) {
        auto* problem = new PlasticityProblem();
        problem->dimension = dim;
        problem->numberOfCones = eta.size();

        // Create M matrix from numpy array
        py::buffer_info M_info = M.request();
        problem->M = NM_create(NM_DENSE, M_info.shape[0], M_info.shape[1]);
        problem->M->matrix0 =
            (double*)malloc(M_info.shape[0] * M_info.shape[1] * sizeof(double));
        std::memcpy(problem->M->matrix0, M_info.ptr,
                    M_info.shape[0] * M_info.shape[1] * sizeof(double));

        // Copy q vector
        py::buffer_info q_info = q.request();
        problem->q = (double*)malloc(q_info.shape[0] * sizeof(double));
        std::memcpy(problem->q, q_info.ptr, q_info.shape[0] * sizeof(double));

        // Create Drucker-Prager model and copy eta/theta
        problem->model_type = PLASTICITY_MODEL_DRUCKER_PRAGER;
        problem->model.drucker_prager =
            (Plasticity_DruckerPrager_model*)malloc(sizeof(Plasticity_DruckerPrager_model));
        py::buffer_info eta_info = eta.request();
        py::buffer_info theta_info = theta.request();
        problem->model.drucker_prager->eta =
            (double*)malloc(eta_info.shape[0] * sizeof(double));
        problem->model.drucker_prager->theta =
            (double*)malloc(theta_info.shape[0] * sizeof(double));
        std::memcpy(problem->model.drucker_prager->eta, eta_info.ptr,
                    eta_info.shape[0] * sizeof(double));
        std::memcpy(problem->model.drucker_prager->theta, theta_info.ptr,
                    theta_info.shape[0] * sizeof(double));

        return problem;
      },
      py::arg("dim"), py::arg("M"), py::arg("q"), py::arg("eta"), py::arg("theta"),
      "Create a PlasticityProblem with Drucker-Prager model",
      py::return_value_policy::take_ownership);

  // Backward compatibility alias
  m.attr("Plasticity2DProblem") = m.attr("PlasticityProblem");

  // Backward compatibility factory function - same as PlasticityProblem_new
  m.def(
      "Plasticity2DProblem_new",
      [](int dim, py::array_t<double> M, py::array_t<double> q, py::array_t<double> eta,
         py::array_t<double> theta) {
        auto* problem = new PlasticityProblem();
        problem->dimension = dim;
        problem->numberOfCones = eta.size();

        // Create M matrix from numpy array
        py::buffer_info M_info = M.request();
        problem->M = NM_create(NM_DENSE, M_info.shape[0], M_info.shape[1]);
        problem->M->matrix0 =
            (double*)malloc(M_info.shape[0] * M_info.shape[1] * sizeof(double));
        std::memcpy(problem->M->matrix0, M_info.ptr,
                    M_info.shape[0] * M_info.shape[1] * sizeof(double));

        // Copy q vector
        py::buffer_info q_info = q.request();
        problem->q = (double*)malloc(q_info.shape[0] * sizeof(double));
        std::memcpy(problem->q, q_info.ptr, q_info.shape[0] * sizeof(double));

        // Create Drucker-Prager model and copy eta/theta
        problem->model_type = PLASTICITY_MODEL_DRUCKER_PRAGER;
        problem->model.drucker_prager =
            (Plasticity_DruckerPrager_model*)malloc(sizeof(Plasticity_DruckerPrager_model));
        py::buffer_info eta_info = eta.request();
        py::buffer_info theta_info = theta.request();
        problem->model.drucker_prager->eta =
            (double*)malloc(eta_info.shape[0] * sizeof(double));
        problem->model.drucker_prager->theta =
            (double*)malloc(theta_info.shape[0] * sizeof(double));
        std::memcpy(problem->model.drucker_prager->eta, eta_info.ptr,
                    eta_info.shape[0] * sizeof(double));
        std::memcpy(problem->model.drucker_prager->theta, theta_info.ptr,
                    theta_info.shape[0] * sizeof(double));

        return problem;
      },
      py::arg("dim"), py::arg("M"), py::arg("q"), py::arg("eta"), py::arg("theta"),
      "Create a PlasticityProblem (backward compatibility)",
      py::return_value_policy::take_ownership);

  // Factory function for creating PlasticityProblem with Von Mises model
  m.def(
      "PlasticityProblem_new_VonMises",
      [](int dim, py::array_t<double> M, py::array_t<double> q, py::array_t<double> sigma_y) {
        auto* problem = new PlasticityProblem();
        problem->dimension = dim;
        problem->numberOfCones = sigma_y.size();
        problem->model_type = PLASTICITY_MODEL_VON_MISES;

        // Create M matrix from numpy array
        py::buffer_info M_info = M.request();
        problem->M = NM_create(NM_DENSE, M_info.shape[0], M_info.shape[1]);
        problem->M->matrix0 =
            (double*)malloc(M_info.shape[0] * M_info.shape[1] * sizeof(double));
        std::memcpy(problem->M->matrix0, M_info.ptr,
                    M_info.shape[0] * M_info.shape[1] * sizeof(double));

        // Copy q vector
        py::buffer_info q_info = q.request();
        problem->q = (double*)malloc(q_info.shape[0] * sizeof(double));
        std::memcpy(problem->q, q_info.ptr, q_info.shape[0] * sizeof(double));

        // Create Von Mises model and copy sigma_y
        problem->model.von_mises =
            (Plasticity_VonMises_model*)malloc(sizeof(Plasticity_VonMises_model));
        py::buffer_info sigma_y_info = sigma_y.request();
        problem->model.von_mises->sigma_y =
            (double*)malloc(sigma_y_info.shape[0] * sizeof(double));
        std::memcpy(problem->model.von_mises->sigma_y, sigma_y_info.ptr,
                    sigma_y_info.shape[0] * sizeof(double));

        return problem;
      },
      py::arg("dim"), py::arg("M"), py::arg("q"), py::arg("sigma_y"),
      "Create a PlasticityProblem with Von Mises model",
      py::return_value_policy::take_ownership);

  // Driver function (forward declaration - actual implementation is in the driver)
  // We need to declare it here to make it available in Python
  // The driver function signature is:
  // int plasticity_2d_driver(PlasticityProblem* problem, double* stress, double*
  // plastic_strain_rate, SolverOptions* options);

  // Display functions
  m.def("plasticity_display", &plasticity_display, "Display plasticity problem",
        py::arg("problem"));

  // Backward compatibility function names
  m.def(
      "mohrCoulomb2D_display", [](PlasticityProblem* problem) { plasticity_display(problem); },
      "Display plasticity problem (backward compatibility)", py::arg("problem"));

  // Note: plasticity_2d_driver is not exposed in Python bindings because
  // it requires numpy array handling for stress and plastic_strain_rate.
  // Users should use the C API or implement a custom wrapper.

  // Solver IDs enum
  py::enum_<PLASTICITY_SOLVER>(solver_ids, "PLASTICITY_SOLVER", "Plasticity solver IDs")
      .value("PLASTICITY_2D_NSGS", PLASTICITY_2D_NSGS, "Non-smooth Gauss-Seidel")
      .value("PLASTICITY_2D_NSGS_GENERIC", PLASTICITY_2D_NSGS_GENERIC,
             "Non-smooth Gauss-Seidel (generic)")
      .value("PLASTICITY_2D_ONECONE_NSN", PLASTICITY_2D_ONECONE_NSN, "Non-smooth Newton")
      .value("PLASTICITY_2D_ONECONE_NSN_GP", PLASTICITY_2D_ONECONE_NSN_GP,
             "Non-smooth Newton GP")
      .value("PLASTICITY_2D_ONECONE_NSN_GP_HYBRID", PLASTICITY_2D_ONECONE_NSN_GP_HYBRID,
             "Hybrid NSN")
      .value("PLASTICITY_2D_ONECONE_ProjectionOnCone", PLASTICITY_2D_ONECONE_ProjectionOnCone,
             "Projection on cone")
      .value("PLASTICITY_2D_ONECONE_ProjectionOnConeWithLocalIteration",
             PLASTICITY_2D_ONECONE_ProjectionOnConeWithLocalIteration,
             "Projection with local iteration")
      .value("PLASTICITY_2D_ONECONE_VONMISES", PLASTICITY_2D_ONECONE_VONMISES,
             "Von Mises radial return")
      .export_values();

  // Also export as module attributes for convenience
  m.attr("PLASTICITY_2D_NSGS") = (int)PLASTICITY_2D_NSGS;
  m.attr("PLASTICITY_2D_NSGS_GENERIC") = (int)PLASTICITY_2D_NSGS_GENERIC;
  m.attr("PLASTICITY_2D_ONECONE_NSN") = (int)PLASTICITY_2D_ONECONE_NSN;
  m.attr("PLASTICITY_2D_ONECONE_NSN_GP") = (int)PLASTICITY_2D_ONECONE_NSN_GP;
  m.attr("PLASTICITY_2D_ONECONE_NSN_GP_HYBRID") = (int)PLASTICITY_2D_ONECONE_NSN_GP_HYBRID;
  m.attr("PLASTICITY_2D_ONECONE_ProjectionOnCone") =
      (int)PLASTICITY_2D_ONECONE_ProjectionOnCone;
  m.attr("PLASTICITY_2D_ONECONE_ProjectionOnConeWithLocalIteration") =
      (int)PLASTICITY_2D_ONECONE_ProjectionOnConeWithLocalIteration;
  m.attr("PLASTICITY_2D_ONECONE_VONMISES") = (int)PLASTICITY_2D_ONECONE_VONMISES;

  // Backward compatibility solver IDs
  m.attr("MOHR_COULOMB_2D_NSGS") = (int)PLASTICITY_2D_NSGS;
  m.attr("SICONOS_MC2D_NSGS") = (int)PLASTICITY_2D_NSGS;
}
