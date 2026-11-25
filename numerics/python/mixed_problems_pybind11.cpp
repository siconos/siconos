/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */

#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>

#include "MixedProblem.h"
#include "NumericsMatrix.h"

namespace py = pybind11;

void wrap_mixed_problems(py::module_& m) {
  // MixedModelType enum
  py::enum_<MixedModelType>(m, "MixedModelType", "Type of mixed problem model")
      .value("UNKNOWN", MIXED_MODEL_UNKNOWN, "Unknown model")
      .value("FRICTION_2D", MIXED_MODEL_FRICTION_2D, "Friction-Contact 2D model")
      .value("FRICTION_3D", MIXED_MODEL_FRICTION_3D, "Friction-Contact 3D model")
      .value("ROLLING_FRICTION_2D", MIXED_MODEL_ROLLING_FRICTION_2D,
             "Rolling Friction-Contact 2D model")
      .value("ROLLING_FRICTION_3D", MIXED_MODEL_ROLLING_FRICTION_3D,
             "Rolling Friction-Contact 3D model")
      .value("VON_MISES", MIXED_MODEL_VON_MISES, "Von Mises plasticity model")
      .value("DRUCKER_PRAGER", MIXED_MODEL_DRUCKER_PRAGER, "Drucker-Prager plasticity model")
      .export_values();

  // Mixed_Friction2D_model
  py::class_<Mixed_Friction2D_model>(
      m, "Mixed_Friction2D_model",
      "Friction-Contact 2D model parameters (friction coefficient)")
      .def(py::init<>())
      .def_readwrite("mu", &Mixed_Friction2D_model::mu,
                     "Friction coefficients (pointer to array)");

  // Mixed_Friction3D_model
  py::class_<Mixed_Friction3D_model>(
      m, "Mixed_Friction3D_model",
      "Friction-Contact 3D model parameters (friction coefficient)")
      .def(py::init<>())
      .def_readwrite("mu", &Mixed_Friction3D_model::mu,
                     "Friction coefficients (pointer to array)");

  // Mixed_RollingFriction2D_model
  py::class_<Mixed_RollingFriction2D_model>(m, "Mixed_RollingFriction2D_model",
                                            "Rolling Friction-Contact 2D model parameters")
      .def(py::init<>())
      .def_readwrite("mu", &Mixed_RollingFriction2D_model::mu,
                     "Friction coefficients (pointer to array)")
      .def_readwrite("mur", &Mixed_RollingFriction2D_model::mur,
                     "Rolling friction coefficients (pointer to array)");

  // Mixed_RollingFriction3D_model
  py::class_<Mixed_RollingFriction3D_model>(m, "Mixed_RollingFriction3D_model",
                                            "Rolling Friction-Contact 3D model parameters")
      .def(py::init<>())
      .def_readwrite("mu", &Mixed_RollingFriction3D_model::mu,
                     "Friction coefficients (pointer to array)")
      .def_readwrite("mur", &Mixed_RollingFriction3D_model::mur,
                     "Rolling friction coefficients (pointer to array)");

  // Mixed_VonMises_model
  py::class_<Mixed_VonMises_model>(m, "Mixed_VonMises_model",
                                   "Von Mises plasticity model parameters (yield stress)")
      .def(py::init<>())
      .def_readwrite("sigma_y", &Mixed_VonMises_model::sigma_y,
                     "Yield stresses (pointer to array)");

  // Mixed_DruckerPrager_model
  py::class_<Mixed_DruckerPrager_model>(m, "Mixed_DruckerPrager_model",
                                        "Drucker-Prager plasticity model parameters")
      .def(py::init<>())
      .def_readwrite("eta", &Mixed_DruckerPrager_model::eta,
                     "Cone coefficients (pointer to array)")
      .def_readwrite("theta", &Mixed_DruckerPrager_model::theta,
                     "Dilatancy coefficients (pointer to array)");

  // MixedProblemBlock
  py::class_<MixedProblemBlock>(m, "MixedProblemBlock", "Block/sub-problem in a mixed problem")
      .def(py::init<>())
      .def_readwrite("start_index", &MixedProblemBlock::start_index,
                     "Start index in global vector")
      .def_readwrite("dimension", &MixedProblemBlock::dimension,
                     "Dimension of this block (2 or 3)")
      .def_readwrite("num_elements", &MixedProblemBlock::num_elements,
                     "Number of elements in this block")
      .def_readwrite("model_type", &MixedProblemBlock::model_type, "Model type enum")
      .def_readwrite("model", &MixedProblemBlock::model, "Model parameters union");

  // MixedProblem
  py::class_<MixedProblem>(
      m, "MixedProblem",
      "Mixed problem structure combining friction-contact and plasticity models")
      .def(py::init<>())
      .def_readwrite("total_dimension", &MixedProblem::total_dimension,
                     "Total dimension of the problem")
      .def_readwrite("num_blocks", &MixedProblem::num_blocks, "Number of blocks/sub-problems")
      .def_readwrite("blocks", &MixedProblem::blocks, "Array of blocks")
      .def_readwrite("M", &MixedProblem::M, "M matrix (shared)")
      .def_readwrite("q", &MixedProblem::q, "q vector (shared)")
      .def(
          "display", [](MixedProblem* problem) { mixedProblem_display(problem); },
          "Display the mixed problem");

  // Factory function for creating an empty MixedProblem
  m.def("MixedProblem_new", &mixedProblem_new, "Create an empty MixedProblem",
        py::return_value_policy::take_ownership);

  // Factory function for creating MixedProblem with data
  m.def(
      "MixedProblem_new_with_data",
      [](int num_blocks, py::array_t<double> M, py::array_t<double> q) {
        py::buffer_info M_info = M.request();
        py::buffer_info q_info = q.request();

        // Create M matrix from numpy array
        NumericsMatrix* M_matrix = NM_create(NM_DENSE, M_info.shape[0], M_info.shape[1]);
        M_matrix->dense = (double*)malloc(M_info.shape[0] * M_info.shape[1] * sizeof(double));
        std::memcpy(M_matrix->dense, M_info.ptr,
                    M_info.shape[0] * M_info.shape[1] * sizeof(double));

        // Copy q vector
        double* q_vec = (double*)malloc(q_info.shape[0] * sizeof(double));
        std::memcpy(q_vec, q_info.ptr, q_info.shape[0] * sizeof(double));

        return mixedProblem_new_with_data(num_blocks, M_matrix, q_vec);
      },
      py::arg("num_blocks"), py::arg("M"), py::arg("q"),
      "Create a MixedProblem with shared M and q", py::return_value_policy::take_ownership);

  // Function to set a block
  m.def("MixedProblem_set_block", &mixedProblem_set_block, "Set block information",
        py::arg("problem"), py::arg("block_index"), py::arg("start_index"),
        py::arg("dimension"), py::arg("num_elements"), py::arg("model_type"));

  // Helper function to set friction 2D params
  m.def(
      "MixedProblem_set_friction2d_params",
      [](MixedProblem* problem, int block_index, py::array_t<double> mu) {
        py::buffer_info mu_info = mu.request();
        double* mu_copy = (double*)malloc(mu_info.shape[0] * sizeof(double));
        std::memcpy(mu_copy, mu_info.ptr, mu_info.shape[0] * sizeof(double));
        return mixedProblem_set_friction2d_params(problem, block_index, mu_copy);
      },
      py::arg("problem"), py::arg("block_index"), py::arg("mu"),
      "Set Friction-Contact 2D parameters for a block");

  // Helper function to set friction 3D params
  m.def(
      "MixedProblem_set_friction3d_params",
      [](MixedProblem* problem, int block_index, py::array_t<double> mu) {
        py::buffer_info mu_info = mu.request();
        double* mu_copy = (double*)malloc(mu_info.shape[0] * sizeof(double));
        std::memcpy(mu_copy, mu_info.ptr, mu_info.shape[0] * sizeof(double));
        return mixedProblem_set_friction3d_params(problem, block_index, mu_copy);
      },
      py::arg("problem"), py::arg("block_index"), py::arg("mu"),
      "Set Friction-Contact 3D parameters for a block");

  // Helper function to set rolling friction 2D params
  m.def(
      "MixedProblem_set_rolling_friction2d_params",
      [](MixedProblem* problem, int block_index, py::array_t<double> mu,
         py::array_t<double> mur) {
        py::buffer_info mu_info = mu.request();
        py::buffer_info mur_info = mur.request();

        double* mu_copy = (double*)malloc(mu_info.shape[0] * sizeof(double));
        double* mur_copy = (double*)malloc(mur_info.shape[0] * sizeof(double));

        std::memcpy(mu_copy, mu_info.ptr, mu_info.shape[0] * sizeof(double));
        std::memcpy(mur_copy, mur_info.ptr, mur_info.shape[0] * sizeof(double));

        return mixedProblem_set_rolling_friction2d_params(problem, block_index, mu_copy,
                                                          mur_copy);
      },
      py::arg("problem"), py::arg("block_index"), py::arg("mu"), py::arg("mur"),
      "Set Rolling Friction-Contact 2D parameters for a block");

  // Helper function to set rolling friction 3D params
  m.def(
      "MixedProblem_set_rolling_friction3d_params",
      [](MixedProblem* problem, int block_index, py::array_t<double> mu,
         py::array_t<double> mur) {
        py::buffer_info mu_info = mu.request();
        py::buffer_info mur_info = mur.request();

        double* mu_copy = (double*)malloc(mu_info.shape[0] * sizeof(double));
        double* mur_copy = (double*)malloc(mur_info.shape[0] * sizeof(double));

        std::memcpy(mu_copy, mu_info.ptr, mu_info.shape[0] * sizeof(double));
        std::memcpy(mur_copy, mur_info.ptr, mur_info.shape[0] * sizeof(double));

        return mixedProblem_set_rolling_friction3d_params(problem, block_index, mu_copy,
                                                          mur_copy);
      },
      py::arg("problem"), py::arg("block_index"), py::arg("mu"), py::arg("mur"),
      "Set Rolling Friction-Contact 3D parameters for a block");

  // Helper function to set Von Mises params
  m.def(
      "MixedProblem_set_vonmises_params",
      [](MixedProblem* problem, int block_index, py::array_t<double> sigma_y) {
        py::buffer_info sigma_y_info = sigma_y.request();
        double* sigma_y_copy = (double*)malloc(sigma_y_info.shape[0] * sizeof(double));
        std::memcpy(sigma_y_copy, sigma_y_info.ptr, sigma_y_info.shape[0] * sizeof(double));
        return mixedProblem_set_vonmises_params(problem, block_index, sigma_y_copy);
      },
      py::arg("problem"), py::arg("block_index"), py::arg("sigma_y"),
      "Set Von Mises parameters for a block");

  // Helper function to set Drucker-Prager params
  m.def(
      "MixedProblem_set_druckerprager_params",
      [](MixedProblem* problem, int block_index, py::array_t<double> eta,
         py::array_t<double> theta) {
        py::buffer_info eta_info = eta.request();
        py::buffer_info theta_info = theta.request();

        double* eta_copy = (double*)malloc(eta_info.shape[0] * sizeof(double));
        double* theta_copy = (double*)malloc(theta_info.shape[0] * sizeof(double));

        std::memcpy(eta_copy, eta_info.ptr, eta_info.shape[0] * sizeof(double));
        std::memcpy(theta_copy, theta_info.ptr, theta_info.shape[0] * sizeof(double));

        return mixedProblem_set_druckerprager_params(problem, block_index, eta_copy,
                                                     theta_copy);
      },
      py::arg("problem"), py::arg("block_index"), py::arg("eta"), py::arg("theta"),
      "Set Drucker-Prager parameters for a block");

  // Utility function to get string representation of model type
  m.def("MixedModelType_to_string", &mixedModelType_to_string,
        "Get string representation of model type", py::arg("model_type"));
}
