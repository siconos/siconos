/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2025 INRIA.
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
#include <pybind11/pytypes.h>
#include <stdint.h>

#include <cstddef>

// #include <cassert>

#include "FrictionContactProblem.h"
#include "FrictionContact_options.h"
#include "GenericMechanicalProblem.h"
#include "GlobalFrictionContactProblem.h"
#include "GlobalRollingFrictionContactProblem.h"
#include "NonSmoothDrivers.h"
#include "NonSmoothSolvers/NonSmoothGaussSeidel_options.h"
#include "NumericsFwd.h"
#include "NumericsMatrix_pybind11.h"
#include "RollingFrictionContactProblem.h"
#include "RollingFrictionContact_options.h"
#include "SolverOptions.h"

namespace py = pybind11;

void display_csc_matrix(const CSparseMatrix* csc) {
  std::cout << "START ... " << csc->nzmax << " " << csc->n << "\n";
  std::cout << "Adress csc->p : " << static_cast<void*>(csc->p) << std::endl;
  std::cout << "Adress  csc->i : " << static_cast<void*>(csc->i) << std::endl;
  std::cout << "Adresd csc->x : " << static_cast<void*>(csc->x) << std::endl;
  std::cout << "p:" << std::endl;
  for (int j = 0; j < csc->n + 1; ++j) {
    std::cout << csc->p[j] << " ";
  }
  std::cout << std::endl;

  std::cout << "i):" << std::endl;
  for (int j = 0; j < csc->nzmax; ++j) {
    std::cout << csc->i[j] << " ";
  }
  std::cout << std::endl;

  std::cout << "(x):" << std::endl;
  for (int j = 0; j < csc->nzmax; ++j) {
    std::cout << csc->x[j] << " ";
  }

  std::cout << std::endl;
}

struct FrictionContactProblemWrapper {
  FrictionContactProblem* problem_{nullptr};
  std::shared_ptr<SparseMatrixMemoryOwner> mem_for_Msparse{nullptr};
  std::shared_ptr<DenseMatrixMemoryOwner> mem_for_Mdense{nullptr};
  py::array_t<double> q_array;
  py::array_t<double> mu_array;

  FrictionContactProblemWrapper(int dimension, int nc, py::object matrix,
                                py::array_t<double> q, py::array_t<double> mu) {
    // Check input vectors and matrices sizes
    if (q.ndim() != 1 || q.shape(0) != dimension * nc)
      throw std::runtime_error("q must be a 1D array of size dim * nc");

    if (mu.ndim() != 1 || mu.shape(0) != nc)
      throw std::runtime_error("mu must be a 1D array of size nc");

    q_array = q;
    mu_array = mu;

    NumericsMatrix* nm;
    if (py::hasattr(matrix, "data") && py::hasattr(matrix, "indices") &&
        py::hasattr(matrix, "indptr")) {
      nm = set_matrix_sparse_new(*this, matrix, mem_for_Msparse);
    } else if (py::isinstance<py::array>(matrix)) {
      nm = set_matrix_dense_new(*this, matrix, mem_for_Mdense);
    } else {
      throw std::runtime_error("Unsupported matrix type");
    }
    problem_ = new FrictionContactProblem;
    problem_->dimension = dimension;
    problem_->numberOfContacts = nc;
    problem_->M = nm;
    problem_->q = static_cast<double*>(q_array.mutable_data());
    problem_->mu = static_cast<double*>(mu_array.mutable_data());
    assert(problem_->M->size0 == problem_->M->size1);
    assert(problem_->M->size0 == dimension * nc);
  }

  ~FrictionContactProblemWrapper() {
    if (problem_) {
      if (problem_->M) {
        if (problem_->M->matrix2) {
          delete problem_->M->matrix2->csc;
          delete problem_->M->matrix2;
        }
        delete problem_->M;
      }
      delete problem_;
    }
  }

  py::object M() const {
    if (problem_->M->storageType == NM_SPARSE) {
      py::module sp = py::module::import("scipy.sparse");
      return sp.attr("csc_matrix")(
          py::make_tuple(mem_for_Msparse->data, mem_for_Msparse->indices,
                         mem_for_Msparse->indptr),
          py::make_tuple(problem_->M->size0, problem_->M->size1));
    } else {
      return mem_for_Mdense->data;
    }
  }

  py::array_t<double> q() const { return q_array; }
  py::array_t<double> mu() const { return mu_array; }

  int solve(py::array_t<double> r, py::array_t<double> w, SolverOptions* solver_options) {
    int info = -1;
    int default_solver_options = 0;
    if (problem_->dimension == 2) {
      if (solver_options == nullptr) {
        solver_options = solver_options_create(SICONOS_FRICTION_2D_NSGS);
        default_solver_options = 1;
      }
      info = fc2d_driver(problem_, r.mutable_data(), w.mutable_data(), solver_options);
      if (default_solver_options) {
        solver_options_delete(solver_options);
        solver_options = nullptr;
      }
    } else if (problem_->dimension == 3) {
      if (solver_options == nullptr) {
        solver_options = solver_options_create(SICONOS_FRICTION_3D_NSGS);
        default_solver_options = 1;
      }
      info = fc3d_driver(problem_, r.mutable_data(), w.mutable_data(), solver_options);
      if (default_solver_options) {
        solver_options_delete(solver_options);
        solver_options = nullptr;
      }
    } else {
      std::cout << "FrictionContactProblem solve: wrong dimension" << std::endl;
    }
    return info;
  };

  void somefunc() {
    // test function ...
    if (problem_->M->storageType == NM_SPARSE) {
      CSparseMatrix* mat = problem_->M->matrix2->csc;
      for (int i = 0; i < mat->nzmax; ++i) {
        mat->x[i] += 10.0;
      }
    } else {
      double* mat = problem_->M->matrix0;
      int size = problem_->M->size0 * problem_->M->size1;
      for (int i = 0; i < size; ++i) {
        mat[i] += 5.0;
      }
    }

    for (int i = 0; i < problem_->dimension * problem_->numberOfContacts; ++i)
      problem_->q[i] += 100.0;
  }
};

struct RollingFrictionContactProblemWrapper {
  RollingFrictionContactProblem* problem{nullptr};
  std::shared_ptr<SparseMatrixMemoryOwner> mem_for_Msparse{nullptr};
  std::shared_ptr<DenseMatrixMemoryOwner> mem_for_Mdense{nullptr};
  py::array_t<double> q_array;
  py::array_t<double> mu_array;
  py::array_t<double> mu_r_array;

  RollingFrictionContactProblemWrapper(int dimension, int nc, py::object matrix,
                                       py::array_t<double> q, py::array_t<double> mu,
                                       py::array_t<double> mu_r) {
    // Check input vectors and matrices sizes
    if (q.ndim() != 1 || q.shape(0) != dimension * nc)
      throw std::runtime_error("q must be a 1D array of size dim * nc");

    if (mu.ndim() != 1 || mu.shape(0) != nc)
      throw std::runtime_error("mu must be a 1D array of size nc");

    if (mu_r.ndim() != 1 || mu_r.shape(0) != nc)
      throw std::runtime_error("mu must be a 1D array of size nc");

    q_array = q;
    mu_array = mu;
    mu_r_array = mu_r;

    NumericsMatrix* nm;
    if (py::hasattr(matrix, "data") && py::hasattr(matrix, "indices") &&
        py::hasattr(matrix, "indptr")) {
      nm = set_matrix_sparse_new(*this, matrix, mem_for_Msparse);
    } else if (py::isinstance<py::array>(matrix)) {
      nm = set_matrix_dense_new(*this, matrix, mem_for_Mdense);
    } else {
      throw std::runtime_error("Unsupported matrix type");
    }

    problem = new RollingFrictionContactProblem;
    problem->dimension = dimension;
    problem->numberOfContacts = nc;
    problem->M = nm;
    problem->q = static_cast<double*>(q_array.mutable_data());
    problem->mu = static_cast<double*>(mu_array.mutable_data());
    problem->mu_r = static_cast<double*>(mu_r_array.mutable_data());
    assert(problem->M->size0 == problem->M->size1);
    assert(problem->M->size0 == dimension * nc);
  }

  ~RollingFrictionContactProblemWrapper() {
    if (problem) {
      if (problem->M) {
        if (problem->M->matrix2) {
          delete problem->M->matrix2->csc;
          delete problem->M->matrix2;
        }
        delete problem->M;
      }
      delete problem;
    }
  }

  py::object M() const {
    if (problem->M->storageType == NM_SPARSE) {
      py::module sp = py::module::import("scipy.sparse");
      return sp.attr("csc_matrix")(
          py::make_tuple(mem_for_Msparse->data, mem_for_Msparse->indices,
                         mem_for_Msparse->indptr),
          py::make_tuple(problem->M->size0, problem->M->size1));
    } else {
      return mem_for_Mdense->data;
    }
  }

  py::array_t<double> q() const { return q_array; }
  py::array_t<double> mu() const { return mu_array; }
  py::array_t<double> mu_r() const { return mu_r_array; }
};

struct GlobalFrictionContactProblemWrapper {
  GlobalFrictionContactProblem* problem{nullptr};
  std::shared_ptr<SparseMatrixMemoryOwner> mem_for_Msparse{nullptr};
  std::shared_ptr<DenseMatrixMemoryOwner> mem_for_Mdense{nullptr};
  std::shared_ptr<SparseMatrixMemoryOwner> mem_for_Hsparse{nullptr};
  std::shared_ptr<DenseMatrixMemoryOwner> mem_for_Hdense{nullptr};
  py::array_t<double> q_array;
  py::array_t<double> mu_array;

  GlobalFrictionContactProblemWrapper(int dimension, int nc, py::object matrixM,
                                      py::object matrixH, py::array_t<double> q,
                                      py::array_t<double> mu) {
    // Check input vectors and matrices sizes
    if (q.ndim() != 1) throw std::runtime_error("q must be a 1D array");

    if (mu.ndim() != 1 || mu.shape(0) != nc)
      throw std::runtime_error("mu must be a 1D array of size nc");

    q_array = q;
    mu_array = mu;

    NumericsMatrix* nm;
    if (py::hasattr(matrixM, "data") && py::hasattr(matrixM, "indices") &&
        py::hasattr(matrixM, "indptr")) {
      nm = set_matrix_sparse_new(*this, matrixM, mem_for_Msparse);
    } else if (py::isinstance<py::array>(matrixM)) {
      nm = set_matrix_dense_new(*this, matrixM, mem_for_Mdense);
    } else {
      throw std::runtime_error("Unsupported matrix type");
    }

    NumericsMatrix* nH;
    if (py::hasattr(matrixH, "data") && py::hasattr(matrixH, "indices") &&
        py::hasattr(matrixH, "indptr")) {
      nH = set_matrix_sparse_new(*this, matrixH, mem_for_Hsparse);
    } else if (py::isinstance<py::array>(matrixH)) {
      nH = set_matrix_dense_new(*this, matrixH, mem_for_Hdense);
    } else {
      throw std::runtime_error("Unsupported matrix type");
    }

    problem = new GlobalFrictionContactProblem;
    problem->dimension = dimension;
    problem->numberOfContacts = nc;
    problem->M = nm;
    problem->H = nH;
    problem->q = static_cast<double*>(q_array.mutable_data());
    problem->mu = static_cast<double*>(mu_array.mutable_data());
    assert(problem->M->size0 == problem->M->size1);
    assert(problem->H->size0 == problem->M->size0);
    assert(problem->H->size1 == dimension * nc);
    assert(q.shape(0) == problem->M->size0);
  }

  GlobalFrictionContactProblemWrapper(std::string fname) {
    problem = globalFrictionContact_new_from_filename(fname.c_str());
  }

  ~GlobalFrictionContactProblemWrapper() {
    if (problem) {
      if (problem->M) {
        if (problem->M->matrix2) {
          delete problem->M->matrix2->csc;
          delete problem->M->matrix2;
        }
        delete problem->M;
      }
      if (problem->H) {
        if (problem->H->matrix2) {
          delete problem->H->matrix2->csc;
          delete problem->H->matrix2;
        }
        delete problem->H;
      }
      delete problem;
    }
  }

  py::object M() const {
    if (problem->M->storageType == NM_SPARSE) {
      py::module sp = py::module::import("scipy.sparse");
      return sp.attr("csc_matrix")(
          py::make_tuple(mem_for_Msparse->data, mem_for_Msparse->indices,
                         mem_for_Msparse->indptr),
          py::make_tuple(problem->M->size0, problem->M->size1));
    } else {
      return mem_for_Mdense->data;
    }
  }
  py::object H() const {
    if (problem->H->storageType == NM_SPARSE) {
      py::module sp = py::module::import("scipy.sparse");
      return sp.attr("csc_matrix")(
          py::make_tuple(mem_for_Hsparse->data, mem_for_Hsparse->indices,
                         mem_for_Hsparse->indptr),
          py::make_tuple(problem->H->size0, problem->H->size1));
    } else {
      return mem_for_Hdense->data;
    }
  }

  py::array_t<double> q() const { return q_array; }
  py::array_t<double> mu() const { return mu_array; }

  void somefunc() {
    // test function ...
    if (problem->M->storageType == NM_SPARSE) {
      CSparseMatrix* mat = problem->M->matrix2->csc;
      for (int i = 0; i < mat->nzmax; ++i) {
        mat->x[i] += 10.0;
      }
    } else {
      double* mat = problem->M->matrix0;
      int size = problem->M->size0 * problem->M->size1;
      for (int i = 0; i < size; ++i) {
        mat[i] += 5.0;
      }
    }

    for (int i = 0; i < problem->dimension * problem->numberOfContacts; ++i)
      problem->q[i] += 100.0;
  }
};

struct GlobalRollingFrictionContactProblemWrapper {
  GlobalRollingFrictionContactProblem* problem{nullptr};
  std::shared_ptr<SparseMatrixMemoryOwner> mem_for_Msparse{nullptr};
  std::shared_ptr<DenseMatrixMemoryOwner> mem_for_Mdense{nullptr};
  std::shared_ptr<SparseMatrixMemoryOwner> mem_for_Hsparse{nullptr};
  std::shared_ptr<DenseMatrixMemoryOwner> mem_for_Hdense{nullptr};
  py::array_t<double> q_array;
  py::array_t<double> mu_array;
  py::array_t<double> mu_r_array;

  GlobalRollingFrictionContactProblemWrapper(int dimension, int nc, py::object matrixM,
                                             py::object matrixH, py::array_t<double> q,
                                             py::array_t<double> mu,
                                             py::array_t<double> mu_r) {
    // Check input vectors and matrices sizes
    if (q.ndim() != 1) throw std::runtime_error("q must be a 1D array");

    if (mu.ndim() != 1 || mu.shape(0) != nc)
      throw std::runtime_error("mu must be a 1D array of size nc");

    if (mu_r.ndim() != 1 || mu_r.shape(0) != nc)
      throw std::runtime_error("mu must be a 1D array of size nc");

    q_array = q;
    mu_array = mu;
    mu_r_array = mu_r;

    NumericsMatrix* nm;
    if (py::hasattr(matrixM, "data") && py::hasattr(matrixM, "indices") &&
        py::hasattr(matrixM, "indptr")) {
      nm = set_matrix_sparse_new(*this, matrixM, mem_for_Msparse);
    } else if (py::isinstance<py::array>(matrixM)) {
      nm = set_matrix_dense_new(*this, matrixM, mem_for_Mdense);
    } else {
      throw std::runtime_error("Unsupported matrix type");
    }

    NumericsMatrix* nH;
    if (py::hasattr(matrixH, "data") && py::hasattr(matrixH, "indices") &&
        py::hasattr(matrixH, "indptr")) {
      nH = set_matrix_sparse_new(*this, matrixH, mem_for_Hsparse);
    } else if (py::isinstance<py::array>(matrixH)) {
      nH = set_matrix_dense_new(*this, matrixH, mem_for_Hdense);
    } else {
      throw std::runtime_error("Unsupported matrix type");
    }

    problem = new GlobalRollingFrictionContactProblem;
    problem->dimension = dimension;
    problem->numberOfContacts = nc;
    problem->M = nm;
    problem->H = nH;
    problem->q = static_cast<double*>(q_array.mutable_data());
    problem->mu = static_cast<double*>(mu_array.mutable_data());
    problem->mu_r = static_cast<double*>(mu_r_array.mutable_data());
    assert(problem->M->size0 == problem->M->size1);
    assert(problem->H->size0 == problem->M->size0);
    assert(problem->H->size1 == dimension * nc);
    assert(q.shape(0) == problem->M->size0);
  }

  ~GlobalRollingFrictionContactProblemWrapper() {
    if (problem) {
      if (problem->M) {
        if (problem->M->matrix2) {
          delete problem->M->matrix2->csc;
          delete problem->M->matrix2;
        }
        delete problem->M;
      }
      if (problem->H) {
        if (problem->H->matrix2) {
          delete problem->H->matrix2->csc;
          delete problem->H->matrix2;
        }
        delete problem->H;
      }
      delete problem;
    }
  }

  py::object M() const {
    if (problem->M->storageType == NM_SPARSE) {
      py::module sp = py::module::import("scipy.sparse");
      return sp.attr("csc_matrix")(
          py::make_tuple(mem_for_Msparse->data, mem_for_Msparse->indices,
                         mem_for_Msparse->indptr),
          py::make_tuple(problem->M->size0, problem->M->size1));
    } else {
      return mem_for_Mdense->data;
    }
  }
  py::object H() const {
    if (problem->H->storageType == NM_SPARSE) {
      py::module sp = py::module::import("scipy.sparse");
      return sp.attr("csc_matrix")(
          py::make_tuple(mem_for_Hsparse->data, mem_for_Hsparse->indices,
                         mem_for_Hsparse->indptr),
          py::make_tuple(problem->H->size0, problem->H->size1));
    } else {
      return mem_for_Hdense->data;
    }
  }

  py::array_t<double> q() const { return q_array; }
  py::array_t<double> mu() const { return mu_array; }
  py::array_t<double> mu_r() const { return mu_r_array; }

  void somefunc() {
    // test function ...
    if (problem->M->storageType == NM_SPARSE) {
      CSparseMatrix* mat = problem->M->matrix2->csc;
      for (int i = 0; i < mat->nzmax; ++i) {
        mat->x[i] += 10.0;
      }
    } else {
      double* mat = problem->M->matrix0;
      int size = problem->M->size0 * problem->M->size1;
      for (int i = 0; i < size; ++i) {
        mat[i] += 5.0;
      }
    }

    for (int i = 0; i < problem->dimension * problem->numberOfContacts; ++i)
      problem->q[i] += 100.0;
  }
};

void wrap_friction_contact(py::module_& m, py::module_& params, py::module_& solver_ids) {
  py::class_<FrictionContactProblemWrapper>(m, "FrictionContactProblem")
      .def(py::init<int, int, py::object, py::array_t<double>, py::array_t<double>>(),
           py::arg("dimension"), py::arg("number_of_contacts"), py::arg("matrix"),
           py::arg("q"), py::arg("mu"))
      .def("M", &FrictionContactProblemWrapper::M)
      .def("q", &FrictionContactProblemWrapper::q)
      .def("mu", &FrictionContactProblemWrapper::mu)
      .def_property_readonly(
          "dimension",
          [](const FrictionContactProblemWrapper& self) { return self.problem_->dimension; })
      .def_property_readonly("numberOfContacts",
                             [](const FrictionContactProblemWrapper& self) {
                               return self.problem_->numberOfContacts;
                             })
      .def("__repr__",
           [](const FrictionContactProblemWrapper& self) {
             std::stringstream ss;
             ss << "<FrictionContactProblemWrapper>\n";
             ss << "  dimension: " << self.problem_->dimension << "\n";
             ss << "  numberOfContacts: " << self.problem_->numberOfContacts << "\n";
             ss << "  M:\n" << py::str(self.M()).cast<std::string>() << "\n";
             ss << "  q: " << py::str(self.q()).cast<std::string>() << "\n";
             ss << "  mu: " << py::str(self.mu()).cast<std::string>() << "\n";
             if (self.problem_->M->matrix2 != nullptr) {
               CSparseMatrix* csc = self.problem_->M->matrix2->csc;
               display_csc_matrix(csc);
             } else {
               ss << "  Dense matrix (matrix0):\n";
               for (int i = 0; i < self.problem_->M->size0; ++i) {
                 for (int j = 0; j < self.problem_->M->size1; ++j) {
                   ss << self.problem_->M->matrix0[i * self.problem_->M->size1 + j] << " ";
                 }
                 ss << "\n";
               }
             }

             return ss.str();
           })
      .def("somefunc", &FrictionContactProblemWrapper::somefunc)
      .def("solve", &FrictionContactProblemWrapper::solve, py::arg("reactions"),
           py::arg("velocities"), py::arg("solver_options") = nullptr);

  py::class_<RollingFrictionContactProblemWrapper>(m, "RollingFrictionContactProblem")
      .def(py::init<int, int, py::object, py::array_t<double>, py::array_t<double>,
                    py::array_t<double>>(),
           py::arg("dimension"), py::arg("number_of_contacts"), py::arg("matrix"),
           py::arg("q"), py::arg("mu"), py::arg("mu_r"))
      .def("M", &RollingFrictionContactProblemWrapper::M)
      .def("q", &RollingFrictionContactProblemWrapper::q)
      .def("mu", &RollingFrictionContactProblemWrapper::mu)
      .def("mu_r", &RollingFrictionContactProblemWrapper::mu_r)
      .def_property_readonly("dimension",
                             [](const RollingFrictionContactProblemWrapper& self) {
                               return self.problem->dimension;
                             })
      .def_property_readonly("numberOfContacts",
                             [](const RollingFrictionContactProblemWrapper& self) {
                               return self.problem->numberOfContacts;
                             })
      .def("__repr__", [](const RollingFrictionContactProblemWrapper& self) {
        std::stringstream ss;
        ss << "<RollingFrictionContactProblemWrapper>\n";
        ss << "  dimension: " << self.problem->dimension << "\n";
        ss << "  numberOfContacts: " << self.problem->numberOfContacts << "\n";
        ss << "  M:\n" << py::str(self.M()).cast<std::string>() << "\n";
        ss << "  q: " << py::str(self.q()).cast<std::string>() << "\n";
        ss << "  mu: " << py::str(self.mu()).cast<std::string>() << "\n";
        if (self.problem->M->matrix2 != nullptr) {
          CSparseMatrix* csc = self.problem->M->matrix2->csc;
          display_csc_matrix(csc);
        } else {
          ss << "  Dense matrix (matrix0):\n";
          for (int i = 0; i < self.problem->M->size0; ++i) {
            for (int j = 0; j < self.problem->M->size1; ++j) {
              ss << self.problem->M->matrix0[i * self.problem->M->size1 + j] << " ";
            }
            ss << "\n";
          }
        }

        return ss.str();
      });

  py::class_<GlobalFrictionContactProblemWrapper>(m, "GlobalFrictionContactProblem")
      .def(py::init<int, int, py::object, py::object, py::array_t<double>,
                    py::array_t<double>>(),
           py::arg("dimension"), py::arg("number_of_contacts"), py::arg("M"), py::arg("H"),
           py::arg("q"), py::arg("mu"), "Build with all operators/parameters (H, M, ...)")
      .def(py::init<std::string>(), py::arg("fname"), "Build from hdf5 file")
      .def("M", &GlobalFrictionContactProblemWrapper::M)
      .def("H", &GlobalFrictionContactProblemWrapper::H)
      .def("q", &GlobalFrictionContactProblemWrapper::q)
      .def("mu", &GlobalFrictionContactProblemWrapper::mu, "Friction coefficient")
      .def_property_readonly("dimension",
                             [](const GlobalFrictionContactProblemWrapper& self) {
                               return self.problem->dimension;
                             })
      .def_property_readonly("numberOfContacts",
                             [](const GlobalFrictionContactProblemWrapper& self) {
                               return self.problem->numberOfContacts;
                             })
      .def("__repr__",
           [](const GlobalFrictionContactProblemWrapper& self) {
             std::stringstream ss;
             ss << "<GlobalFrictionContactProblemWrapper>\n";
             ss << "  dimension: " << self.problem->dimension << "\n";
             ss << "  numberOfContacts: " << self.problem->numberOfContacts << "\n";
             ss << "  M:\n" << py::str(self.M()).cast<std::string>() << "\n";
             ss << "  H:\n" << py::str(self.H()).cast<std::string>() << "\n";
             ss << "  q: " << py::str(self.q()).cast<std::string>() << "\n";
             ss << "  mu: " << py::str(self.mu()).cast<std::string>() << "\n";
             if (self.problem->M->matrix2 != nullptr) {
               CSparseMatrix* csc = self.problem->M->matrix2->csc;
               display_csc_matrix(csc);
             } else {
               ss << "  Dense matrix (matrix0):\n";
               for (int i = 0; i < self.problem->M->size0; ++i) {
                 for (int j = 0; j < self.problem->M->size1; ++j) {
                   ss << self.problem->M->matrix0[i * self.problem->M->size1 + j] << " ";
                 }
                 ss << "\n";
               }
             }

             return ss.str();
           })
      .def("somefunc", &GlobalFrictionContactProblemWrapper::somefunc);

  py::class_<GlobalRollingFrictionContactProblemWrapper>(m,
                                                         "GlobalRollingFrictionContactProblem")
      .def(py::init<int, int, py::object, py::object, py::array_t<double>, py::array_t<double>,
                    py::array_t<double>>(),
           py::arg("dimension"), py::arg("number_of_contacts"), py::arg("M"), py::arg("H"),
           py::arg("q"), py::arg("mu"), py::arg("mu_r"))
      .def("M", &GlobalRollingFrictionContactProblemWrapper::M)
      .def("H", &GlobalRollingFrictionContactProblemWrapper::H)
      .def("q", &GlobalRollingFrictionContactProblemWrapper::q)
      .def("mu", &GlobalRollingFrictionContactProblemWrapper::mu)
      .def("mu_r", &GlobalRollingFrictionContactProblemWrapper::mu_r)
      .def_property_readonly("dimension",
                             [](const GlobalRollingFrictionContactProblemWrapper& self) {
                               return self.problem->dimension;
                             })
      .def_property_readonly("numberOfContacts",
                             [](const GlobalRollingFrictionContactProblemWrapper& self) {
                               return self.problem->numberOfContacts;
                             })
      .def("__repr__",
           [](const GlobalRollingFrictionContactProblemWrapper& self) {
             std::stringstream ss;
             ss << "<GlobalFrictionContactProblemWrapper>\n";
             ss << "  dimension: " << self.problem->dimension << "\n";
             ss << "  numberOfContacts: " << self.problem->numberOfContacts << "\n";
             ss << "  M:\n" << py::str(self.M()).cast<std::string>() << "\n";
             ss << "  H:\n" << py::str(self.H()).cast<std::string>() << "\n";
             ss << "  q: " << py::str(self.q()).cast<std::string>() << "\n";
             ss << "  mu: " << py::str(self.mu()).cast<std::string>() << "\n";
             ss << "  mu_r: " << py::str(self.mu_r()).cast<std::string>() << "\n";
             if (self.problem->M->matrix2 != nullptr) {
               CSparseMatrix* csc = self.problem->M->matrix2->csc;
               display_csc_matrix(csc);
             } else {
               ss << "  Dense matrix (matrix0):\n";
               for (int i = 0; i < self.problem->M->size0; ++i) {
                 for (int j = 0; j < self.problem->M->size1; ++j) {
                   ss << self.problem->M->matrix0[i * self.problem->M->size1 + j] << " ";
                 }
                 ss << "\n";
               }
             }

             return ss.str();
           })
      .def("somefunc", &GlobalRollingFrictionContactProblemWrapper::somefunc);

  // Exposing the FRICTION_SOLVER enum using .value
  py::enum_<FRICTION_SOLVER>(solver_ids, "FRICTION_SOLVER",
                             "Friction solver types for 2D and 3D frictional contact problems")
      // 2D Frictional Contact solvers
      .value("SICONOS_FRICTION_2D_NSGS", FRICTION_SOLVER::SICONOS_FRICTION_2D_NSGS,
             "2D Non-smooth Gauss Seidel")
      .value("SICONOS_FRICTION_2D_CPG", FRICTION_SOLVER::SICONOS_FRICTION_2D_CPG,
             "2D CPG solver")
      .value("SICONOS_FRICTION_2D_LEMKE", FRICTION_SOLVER::SICONOS_FRICTION_2D_LEMKE,
             "2D Lemke solver")
      .value("SICONOS_FRICTION_2D_ENUM", FRICTION_SOLVER::SICONOS_FRICTION_2D_ENUM,
             "2D Friction solver enum")

      // 3D frictional contact solvers on local formulation
      .value("SICONOS_FRICTION_3D_NSGS", FRICTION_SOLVER::SICONOS_FRICTION_3D_NSGS,
             "3D Non-smooth Gauss Seidel")
      .value("SICONOS_FRICTION_3D_NSGSV", FRICTION_SOLVER::SICONOS_FRICTION_3D_NSGSV,
             "3D Non-smooth Gauss Seidel-velocity")
      .value("SICONOS_FRICTION_3D_PROX", FRICTION_SOLVER::SICONOS_FRICTION_3D_PROX,
             "3D Proximal solver")
      .value("SICONOS_FRICTION_3D_TFP", FRICTION_SOLVER::SICONOS_FRICTION_3D_TFP,
             "3D Tresca solver")
      .value("SICONOS_FRICTION_3D_NSN_AC", FRICTION_SOLVER::SICONOS_FRICTION_3D_NSN_AC,
             "3D Non-smooth Newton Alart-Curnier")
      .value("SICONOS_FRICTION_3D_DSFP", FRICTION_SOLVER::SICONOS_FRICTION_3D_DSFP,
             "3D De Saxce fixed point")
      .value("SICONOS_FRICTION_3D_VI_FPP", FRICTION_SOLVER::SICONOS_FRICTION_3D_VI_FPP,
             "3D VI formulation, fixed point projection")
      .value("SICONOS_FRICTION_3D_VI_EG", FRICTION_SOLVER::SICONOS_FRICTION_3D_VI_EG,
             "3D VI formulation, Extra-gradient")
      .value("SICONOS_FRICTION_3D_HP", FRICTION_SOLVER::SICONOS_FRICTION_3D_HP,
             "3D Hyperplane projection")
      .value("SICONOS_FRICTION_3D_NCPGlockerFBFixedPoint",
             FRICTION_SOLVER::SICONOS_FRICTION_3D_NCPGlockerFBFixedPoint,
             "3D Fischer Burmeister fixed point")
      .value("SICONOS_FRICTION_3D_FPP", FRICTION_SOLVER::SICONOS_FRICTION_3D_FPP,
             "3D Fixed point projection")
      .value("SICONOS_FRICTION_3D_EG", FRICTION_SOLVER::SICONOS_FRICTION_3D_EG,
             "3D Extra-gradient")
      .value("SICONOS_FRICTION_3D_NSN_FB", FRICTION_SOLVER::SICONOS_FRICTION_3D_NSN_FB,
             "3D Non-smooth Newton Fischer Burmeister")
      .value("SICONOS_FRICTION_3D_GAMS_PATH", FRICTION_SOLVER::SICONOS_FRICTION_3D_GAMS_PATH,
             "3D GAMS/Path (Ferris)")
      .value("SICONOS_FRICTION_3D_GAMS_PATHVI",
             FRICTION_SOLVER::SICONOS_FRICTION_3D_GAMS_PATHVI, "3D GAMS/Path VI formulation")

      // 3D Frictional Contact solvers for one contact (used mainly inside NSGS solvers)
      .value("SICONOS_ONECONE_NSN", FRICTION_SOLVER::SICONOS_ONECONE_NSN,
             "3D One contact Non-smooth Newton Alart-Curnier")
      .value("SICONOS_ONECONE_NSN_GP", FRICTION_SOLVER::SICONOS_ONECONE_NSN_GP,
             "3D One contact Non-smooth Newton Alart-Curnier, damped")
      .value("SICONOS_ONECONE_ProjectionOnCone",
             FRICTION_SOLVER::SICONOS_ONECONE_ProjectionOnCone,
             "3D One contact Projection on cone")
      .value("SICONOS_ONECONE_ProjectionOnConeWithLocalIteration",
             FRICTION_SOLVER::SICONOS_ONECONE_ProjectionOnConeWithLocalIteration,
             "3D One contact Projection on cone with local iteration")
      .value("SICONOS_ONECONE_ProjectionOnConeWithRegularization",
             FRICTION_SOLVER::SICONOS_ONECONE_ProjectionOnConeWithRegularization,
             "3D One contact Projection on cone with regularization")
      .value("SICONOS_ONECONE_ProjectionOnConeWithDiagonalization",
             FRICTION_SOLVER::SICONOS_ONECONE_ProjectionOnConeWithDiagonalization,
             "3D One contact Projection on cone with diagonalization")
      .value("SICONOS_ONECONE_ProjectionOnCone_velocity",
             FRICTION_SOLVER::SICONOS_ONECONE_ProjectionOnCone_velocity,
             "3D One contact Projection on cone, velocity")
      .value("SICONOS_ONECONE_NSN_GP_HYBRID", FRICTION_SOLVER::SICONOS_ONECONE_NSN_GP_HYBRID,
             "3D Frictional contact hybrid solver for one contact")
      .value("SICONOS_FRICTION_3D_VI_FPP_Cylinder",
             FRICTION_SOLVER::SICONOS_FRICTION_3D_VI_FPP_Cylinder,
             "3D VI formulation, Fixed Point Projection Cylinder")
      .value("SICONOS_FRICTION_3D_CONVEXQP_PG_CYLINDER",
             FRICTION_SOLVER::SICONOS_FRICTION_3D_CONVEXQP_PG_CYLINDER,
             "3D Convex QP Projection on Cylinder")
      .value("SICONOS_FRICTION_3D_ACLMFP", FRICTION_SOLVER::SICONOS_FRICTION_3D_ACLMFP,
             "Alart-Curnier fixed point, local formulation")
      .value("SICONOS_FRICTION_3D_SOCLCP", FRICTION_SOLVER::SICONOS_FRICTION_3D_SOCLCP,
             "Second-order Cone LCP, local formulation")
      .value("SICONOS_FRICTION_3D_GAMS_LCP_PATH",
             FRICTION_SOLVER::SICONOS_FRICTION_3D_GAMS_LCP_PATH,
             "GAMS/PATH (Ferris) LCP, local formulation")
      .value("SICONOS_FRICTION_3D_GAMS_LCP_PATHVI",
             FRICTION_SOLVER::SICONOS_FRICTION_3D_GAMS_LCP_PATHVI,
             "VI formulation, GAMS/PATH (Ferris) LCP, local formulation")
      .value("SICONOS_FRICTION_3D_NSN_NM", FRICTION_SOLVER::SICONOS_FRICTION_3D_NSN_NM,
             "Non-smooth Newton, natural map, local formulation")
      .value("SICONOS_FRICTION_3D_NSN_AC_NEW", FRICTION_SOLVER::SICONOS_FRICTION_3D_NSN_AC_NEW,
             "NSN AC Test solver")
      .value("SICONOS_FRICTION_3D_PFP", FRICTION_SOLVER::SICONOS_FRICTION_3D_PFP,
             "Panagiotopoulos, fixed point, local formulation")
      .value("SICONOS_FRICTION_3D_ADMM", FRICTION_SOLVER::SICONOS_FRICTION_3D_ADMM,
             "ADMM local formulation")

      // Fischer Burmeister/Path, Glocker formulation, one contact solver
      .value("SICONOS_FRICTION_3D_NCPGlockerFBPATH",
             FRICTION_SOLVER::SICONOS_FRICTION_3D_NCPGlockerFBPATH,
             "3D Fischer Burmeister, Glocker formulation")
      .value("SICONOS_FRICTION_3D_NCPGlockerFBNewton",
             FRICTION_SOLVER::SICONOS_FRICTION_3D_NCPGlockerFBNewton,
             "3D Fischer Burmeister, Glocker formulation, Newton")
      .value("SICONOS_FRICTION_3D_ONECONTACT_QUARTIC",
             FRICTION_SOLVER::SICONOS_FRICTION_3D_ONECONTACT_QUARTIC,
             "3D One contact Quartic solver")
      .value("SICONOS_FRICTION_3D_ONECONTACT_QUARTIC_NU",
             FRICTION_SOLVER::SICONOS_FRICTION_3D_ONECONTACT_QUARTIC_NU,
             "3D One contact Quartic solver, NU")
      .value("SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCylinder",
             FRICTION_SOLVER::SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCylinder,
             "3D One contact Projection on cylinder")
      .value("SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCylinderWithLocalIteration",
             FRICTION_SOLVER::
                 SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCylinderWithLocalIteration,
             "3D One contact Projection on cylinder with local iteration")

      /** 3D Frictional contact local solvers on global formulation */
      .value("SICONOS_GLOBAL_FRICTION_3D_NSGS_WR",
             FRICTION_SOLVER::SICONOS_GLOBAL_FRICTION_3D_NSGS_WR,
             "Global 3D Frictional Contact solver NSGS (wrapped)")
      .value("SICONOS_GLOBAL_FRICTION_3D_NSGSV_WR",
             FRICTION_SOLVER::SICONOS_GLOBAL_FRICTION_3D_NSGSV_WR,
             "Global 3D Frictional Contact solver NSGSV (wrapped)")
      .value("SICONOS_GLOBAL_FRICTION_3D_PROX_WR",
             FRICTION_SOLVER::SICONOS_GLOBAL_FRICTION_3D_PROX_WR,
             "Global 3D Frictional Contact solver PROX (wrapped)")
      .value("SICONOS_GLOBAL_FRICTION_3D_DSFP_WR",
             FRICTION_SOLVER::SICONOS_GLOBAL_FRICTION_3D_DSFP_WR,
             "Global 3D Frictional Contact solver DSFP (wrapped)")
      .value("SICONOS_GLOBAL_FRICTION_3D_TFP_WR",
             FRICTION_SOLVER::SICONOS_GLOBAL_FRICTION_3D_TFP_WR,
             "Global 3D Frictional contact solver TFP (wrapped)")
      .value("SICONOS_GLOBAL_FRICTION_3D_NSGS",
             FRICTION_SOLVER::SICONOS_GLOBAL_FRICTION_3D_NSGS,
             "Global 3D Frictional contact solver NSGS")
      .value("SICONOS_GLOBAL_FRICTION_3D_NSN_AC_WR",
             FRICTION_SOLVER::SICONOS_GLOBAL_FRICTION_3D_NSN_AC_WR,
             "Global 3D Non-smooth Newton Alart-Curnier solver (wrapped)")
      .value("SICONOS_GLOBAL_FRICTION_3D_NSN_AC",
             FRICTION_SOLVER::SICONOS_GLOBAL_FRICTION_3D_NSN_AC,
             "Global 3D Non-smooth Newton Alart-Curnier solver")
      .value("SICONOS_GLOBAL_FRICTION_3D_GAMS_PATH",
             FRICTION_SOLVER::SICONOS_GLOBAL_FRICTION_3D_GAMS_PATH,
             "Global 3D GAMS/Path (Ferris) solver")
      .value("SICONOS_GLOBAL_FRICTION_3D_GAMS_PATHVI",
             FRICTION_SOLVER::SICONOS_GLOBAL_FRICTION_3D_GAMS_PATHVI,
             "Global 3D GAMS/Path VI solver")

      /** VI formulations */
      .value("SICONOS_GLOBAL_FRICTION_3D_VI_FPP",
             FRICTION_SOLVER::SICONOS_GLOBAL_FRICTION_3D_VI_FPP,
             "Global 3D VI formulation, Fixed Point Projection")
      .value("SICONOS_GLOBAL_FRICTION_3D_VI_EG",
             FRICTION_SOLVER::SICONOS_GLOBAL_FRICTION_3D_VI_EG,
             "Global 3D VI formulation, Extra-gradient")
      .value("SICONOS_GLOBAL_FRICTION_3D_ACLMFP",
             FRICTION_SOLVER::SICONOS_GLOBAL_FRICTION_3D_ACLMFP,
             "Global 3D Alart-Curnier fixed point")
      .value("SICONOS_GLOBAL_FRICTION_3D_ADMM",
             FRICTION_SOLVER::SICONOS_GLOBAL_FRICTION_3D_ADMM, "Global 3D ADMM solver")
      .value("SICONOS_GLOBAL_FRICTION_3D_ADMM_WR",
             FRICTION_SOLVER::SICONOS_GLOBAL_FRICTION_3D_ADMM_WR,
             "Global 3D ADMM solver (wrapped)")
      .value("SICONOS_GLOBAL_FRICTION_3D_IPM", FRICTION_SOLVER::SICONOS_GLOBAL_FRICTION_3D_IPM,
             "Global 3D Interior-point method solver")
      .export_values();

  // Rolling Friction solver enum (separate from FRICTION_SOLVER)
  py::enum_<ROLLING_FRICTION_SOLVER>(solver_ids, "ROLLING_FRICTION_SOLVER_enum",
                                     "Rolling friction solver IDs")
      // Rolling Friction solvers (3D and 2D)
      .value("SICONOS_ROLLING_FRICTION_3D_NSGS",
             ROLLING_FRICTION_SOLVER::SICONOS_ROLLING_FRICTION_3D_NSGS,
             "3D Non-smooth Gauss Seidel, local formulation")
      .value("SICONOS_ROLLING_FRICTION_3D_ADMM",
             ROLLING_FRICTION_SOLVER::SICONOS_ROLLING_FRICTION_3D_ADMM,
             "3D Rolling friction ADMM solver")

      // Rolling friction solvers for 2D problems
      .value("SICONOS_ROLLING_FRICTION_2D_NSGS",
             ROLLING_FRICTION_SOLVER::SICONOS_ROLLING_FRICTION_2D_NSGS,
             "2D Non-smooth Gauss Seidel, local formulation")
      .value("SICONOS_ROLLING_FRICTION_2D_ONECONTACT_ProjectionOnCone",
             ROLLING_FRICTION_SOLVER::SICONOS_ROLLING_FRICTION_2D_ONECONTACT_ProjectionOnCone,
             "2D Rolling friction one contact Projection on Cone")
      .value("SICONOS_ROLLING_FRICTION_2D_ONECONTACT_ProjectionOnConeWithLocalIteration",
             ROLLING_FRICTION_SOLVER::
                 SICONOS_ROLLING_FRICTION_2D_ONECONTACT_ProjectionOnConeWithLocalIteration,
             "2D Rolling friction one contact Projection on Cone with local iteration")

      // Global Rolling Friction solvers for 3D problems
      .value("SICONOS_GLOBAL_ROLLING_FRICTION_3D_NSGS_WR",
             ROLLING_FRICTION_SOLVER::SICONOS_GLOBAL_ROLLING_FRICTION_3D_NSGS_WR,
             "Global 3D Rolling friction solver NSGS (wrapped)")
      .value("SICONOS_GLOBAL_ROLLING_FRICTION_3D_IPM",
             ROLLING_FRICTION_SOLVER::SICONOS_GLOBAL_ROLLING_FRICTION_3D_IPM,
             "Global 3D Rolling friction solver IPM")
      .export_values();

  py::enum_<SICONOS_FRICTION_3D_IPARAM>(params, "SICONOS_FRICTION_3D_IPARAM_enum")
      .value("SICONOS_FRICTION_3D_IPARAM_INTERNAL_ERROR_STRATEGY",
             SICONOS_FRICTION_3D_IPARAM_INTERNAL_ERROR_STRATEGY)
      .value("SICONOS_FRICTION_3D_IPARAM_RESCALING", SICONOS_FRICTION_3D_IPARAM_RESCALING)
      .value("SICONOS_FRICTION_3D_IPARAM_RESCALING_CONE",
             SICONOS_FRICTION_3D_IPARAM_RESCALING_CONE)
      .value("SICONOS_FRICTION_3D_CURRENT_CONTACT_NUMBER",
             SICONOS_FRICTION_3D_CURRENT_CONTACT_NUMBER)
      .value("SICONOS_FRICTION_3D_NUMBER_OF_CONTACTS", SICONOS_FRICTION_3D_NUMBER_OF_CONTACTS)
      .export_values();

  py::enum_<SICONOS_NSGS_INTERNAL_ERROR_STRATEGY>(params,
                                                  "SICONOS_NSGS_INTERNAL_ERROR_STRATEGY_enum")
      .value("SICONOS_NSGS_INTERNAL_ERROR_STRATEGY_ADAPTIVE",
             SICONOS_NSGS_INTERNAL_ERROR_STRATEGY_ADAPTIVE)
      .value("SICONOS_NSGS_INTERNAL_ERROR_STRATEGY_GIVEN_VALUE",
             SICONOS_NSGS_INTERNAL_ERROR_STRATEGY_GIVEN_VALUE)
      .value("SICONOS_NSGS_INTERNAL_ERROR_STRATEGY_ADAPTIVE_N_CONTACT",
             SICONOS_NSGS_INTERNAL_ERROR_STRATEGY_ADAPTIVE_N_CONTACT)
      .export_values();

  py::enum_<SICONOS_FRICTION_3D_RESCALING>(params, "SICONOS_FRICTION_3D_RESCALING_enum")
      .value("SICONOS_FRICTION_3D_RESCALING_NO", SICONOS_FRICTION_3D_RESCALING_NO)
      .value("SICONOS_FRICTION_3D_RESCALING_SCALAR", SICONOS_FRICTION_3D_RESCALING_SCALAR)
      .value("SICONOS_FRICTION_3D_RESCALING_BALANCING_M",
             SICONOS_FRICTION_3D_RESCALING_BALANCING_M)
      .value("SICONOS_FRICTION_3D_RESCALING_BALANCING_MH",
             SICONOS_FRICTION_3D_RESCALING_BALANCING_MH)
      .value("SICONOS_FRICTION_3D_RESCALING_BALANCING_MHHT",
             SICONOS_FRICTION_3D_RESCALING_BALANCING_MHHT)
      .export_values();

  py::enum_<SICONOS_FRICTION_3D_RESCALING_CONE>(params,
                                                "SICONOS_FRICTION_3D_RESCALING_CONE_enum")
      .value("SICONOS_FRICTION_3D_RESCALING_CONE_NO", SICONOS_FRICTION_3D_RESCALING_CONE_NO)
      .value("SICONOS_FRICTION_3D_RESCALING_CONE_YES", SICONOS_FRICTION_3D_RESCALING_CONE_YES)
      .export_values();


  py::enum_<SICONOS_NSGS_IPARAM>(params, "SICONOS_NSGS_IPARAM_enum")
      .value("SICONOS_NSGS_RELAXATION", SICONOS_NSGS_RELAXATION)
      .value("SICONOS_NSGS_SHUFFLE", SICONOS_NSGS_SHUFFLE)
      .value("SICONOS_NSGS_SHUFFLE_SEED", SICONOS_NSGS_SHUFFLE_SEED)
      .value("SICONOS_NSGS_FREEZING_CONTACT", SICONOS_NSGS_FREEZING_CONTACT)
      .value("SICONOS_NSGS_FILTER_LOCAL_SOLUTION", SICONOS_NSGS_FILTER_LOCAL_SOLUTION)
      .value("SICONOS_NSGS_ERROR_EVALUATION_TYPE", SICONOS_NSGS_ERROR_EVALUATION_TYPE)
      .value("SICONOS_NSGS_ERROR_EVALUATION_FREQUENCY",
             SICONOS_NSGS_ERROR_EVALUATION_FREQUENCY)
      .export_values();

  py::enum_<SICONOS_NSGS_DPARAM>(params, "SICONOS_NSGS_DPARAM_enum")
      .value("SICONOS_NSGS_RELAXATION_VALUE", SICONOS_NSGS_RELAXATION_VALUE)
      .value("SICONOS_NSGS_INTERNAL_ERROR_RATIO", SICONOS_NSGS_INTERNAL_ERROR_RATIO)
      .export_values();

  py::enum_<enum SICONOS_NSGS_LOCALSOLVER_USE_TRIVIAL_SOLUTION>(
      params, "SICONOS_NSGS_LOCALSOLVER_USE_TRIVIAL_SOLUTION_enum")
      .value("SICONOS_NSGS_LOCALSOLVER_USE_TRIVIAL_SOLUTION_FALSE",
             SICONOS_NSGS_LOCALSOLVER_USE_TRIVIAL_SOLUTION_FALSE)
      .value("SICONOS_NSGS_LOCALSOLVER_USE_TRIVIAL_SOLUTION_TRUE",
             SICONOS_NSGS_LOCALSOLVER_USE_TRIVIAL_SOLUTION_TRUE)
      .export_values();

  py::enum_<SICONOS_NSGS_ERROR_EVALUATION>(params, "SICONOS_NSGS_ERROR_EVALUATION_enum")
      .value("SICONOS_NSGS_ERROR_EVALUATION_FULL", SICONOS_NSGS_ERROR_EVALUATION_FULL)
      .value("SICONOS_NSGS_ERROR_EVALUATION_LIGHT", SICONOS_NSGS_ERROR_EVALUATION_LIGHT)
      .value("SICONOS_NSGS_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL",
             SICONOS_NSGS_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL)
      .value("SICONOS_NSGS_ERROR_EVALUATION_ADAPTIVE", SICONOS_NSGS_ERROR_EVALUATION_ADAPTIVE)
      .export_values();

  py::enum_<enum SICONOS_NSGS_SHUFFLE>(params, "SICONOS_NSGS_SHUFFLE_enum")
      .value("SICONOS_NSGS_SHUFFLE_FALSE", SICONOS_NSGS_SHUFFLE_FALSE)
      .value("SICONOS_NSGS_SHUFFLE_TRUE", SICONOS_NSGS_SHUFFLE_TRUE)
      .value("SICONOS_NSGS_SHUFFLE_EACH_LOOP", SICONOS_NSGS_SHUFFLE_EACH_LOOP)
      .export_values();

  py::enum_<enum SICONOS_NSGS_RELAXATION>(params, "SICONOS_NSGS_RELAXATION_enum")
      .value("SICONOS_NSGS_RELAXATION_FALSE", SICONOS_NSGS_RELAXATION_FALSE)
      .value("SICONOS_NSGS_RELAXATION_TRUE", SICONOS_NSGS_RELAXATION_TRUE)
      .export_values();

  py::enum_<enum SICONOS_NSGS_FILTER_LOCAL_SOLUTION>(params,
                                                     "SICONOS_NSGS_FILTER_LOCAL_SOLUTION_enum")
      .value("SICONOS_NSGS_FILTER_LOCAL_SOLUTION_FALSE",
             SICONOS_NSGS_FILTER_LOCAL_SOLUTION_FALSE)
      .value("SICONOS_NSGS_FILTER_LOCAL_SOLUTION_TRUE",
             SICONOS_NSGS_FILTER_LOCAL_SOLUTION_TRUE)
      .export_values();

  py::enum_<SICONOS_FRICTION_3D_NSN_IPARAM>(params, "SICONOS_FRICTION_3D_NSN_IPARAM_enum")
      .value("SICONOS_FRICTION_3D_NSN_RHO_STRATEGY", SICONOS_FRICTION_3D_NSN_RHO_STRATEGY)
      .value("SICONOS_FRICTION_3D_NSN_FORMULATION", SICONOS_FRICTION_3D_NSN_FORMULATION)
      .value("SICONOS_FRICTION_3D_NSN_LINESEARCH", SICONOS_FRICTION_3D_NSN_LINESEARCH)
      .value("SICONOS_FRICTION_3D_NSN_LINESEARCH_MAX_ITER",
             SICONOS_FRICTION_3D_NSN_LINESEARCH_MAX_ITER)
      .value("SICONOS_FRICTION_3D_NSN_LINEAR_SOLVER", SICONOS_FRICTION_3D_NSN_LINEAR_SOLVER)
      .value("SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY",
             SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY)
      .value("SICONOS_FRICTION_3D_NSN_HYBRID_MAX_LOOP",
             SICONOS_FRICTION_3D_NSN_HYBRID_MAX_LOOP)
      .value("SICONOS_FRICTION_3D_NSN_HYBRID_MAX_ITER",
             SICONOS_FRICTION_3D_NSN_HYBRID_MAX_ITER)
      .value("SICONOS_FRICTION_3D_NSN_MEMORY_ALLOCATED",
             SICONOS_FRICTION_3D_NSN_MEMORY_ALLOCATED)
      .value("SICONOS_FRICTION_3D_NSN_MPI_COM", SICONOS_FRICTION_3D_NSN_MPI_COM)
      .export_values();

  py::enum_<SICONOS_FC3D_NSN_LINEAR_SOLVER>(params, "SICONOS_FC3D_NSN_LINEAR_SOLVER_enum")
      .value("SICONOS_FRICTION_3D_NSN_USE_CSLUSOL", SICONOS_FRICTION_3D_NSN_USE_CSLUSOL)
      .value("SICONOS_FRICTION_3D_NSN_USE_MUMPS", SICONOS_FRICTION_3D_NSN_USE_MUMPS)
      .export_values();

  py::enum_<SICONOS_FRICTION_3D_NSN_DPARAM>(params, "SICONOS_FRICTION_3D_NSN_DPARAM_enum")
      .value("SICONOS_FRICTION_3D_NSN_RHO", SICONOS_FRICTION_3D_NSN_RHO)
      .export_values();

  py::enum_<enum SICONOS_FRICTION_3D_NSN_RHO_STRATEGY>(
      params, "SICONOS_FRICTION_3D_NSN_RHO_STRATEGY_enum")
      .value("SICONOS_FRICTION_3D_NSN_FORMULATION_RHO_STRATEGY_CONSTANT",
             SICONOS_FRICTION_3D_NSN_FORMULATION_RHO_STRATEGY_CONSTANT)
      .value("SICONOS_FRICTION_3D_NSN_FORMULATION_RHO_STRATEGY_SPECTRAL_NORM",
             SICONOS_FRICTION_3D_NSN_FORMULATION_RHO_STRATEGY_SPECTRAL_NORM)
      .value("SICONOS_FRICTION_3D_NSN_FORMULATION_RHO_STRATEGY_SPLIT_SPECTRAL_NORM_COND",
             SICONOS_FRICTION_3D_NSN_FORMULATION_RHO_STRATEGY_SPLIT_SPECTRAL_NORM_COND)
      .value("SICONOS_FRICTION_3D_NSN_FORMULATION_RHO_STRATEGY_SPLIT_SPECTRAL_NORM",
             SICONOS_FRICTION_3D_NSN_FORMULATION_RHO_STRATEGY_SPLIT_SPECTRAL_NORM)
      .value("SICONOS_FRICTION_3D_NSN_FORMULATION_RHO_STRATEGY_ADAPTIVE",
             SICONOS_FRICTION_3D_NSN_FORMULATION_RHO_STRATEGY_ADAPTIVE)
      .export_values();

  py::enum_<enum SICONOS_FRICTION_3D_NSN_FORMULATION>(
      params, "SICONOS_FRICTION_3D_NSN_FORMULATION_enum")
      .value("SICONOS_FRICTION_3D_NSN_FORMULATION_ALARTCURNIER_STD",
             SICONOS_FRICTION_3D_NSN_FORMULATION_ALARTCURNIER_STD)
      .value("SICONOS_FRICTION_3D_NSN_FORMULATION_JEANMOREAU_STD",
             SICONOS_FRICTION_3D_NSN_FORMULATION_JEANMOREAU_STD)
      .value("SICONOS_FRICTION_3D_NSN_FORMULATION_ALARTCURNIER_GENERATED",
             SICONOS_FRICTION_3D_NSN_FORMULATION_ALARTCURNIER_GENERATED)
      .value("SICONOS_FRICTION_3D_NSN_FORMULATION_JEANMOREAU_GENERATED",
             SICONOS_FRICTION_3D_NSN_FORMULATION_JEANMOREAU_GENERATED)
      .value("SICONOS_FRICTION_3D_NSN_FORMULATION_NULL",
             SICONOS_FRICTION_3D_NSN_FORMULATION_NULL)
      .export_values();

  py::enum_<enum SICONOS_FRICTION_3D_NSN_LINESEARCH>(params,
                                                     "SICONOS_FRICTION_3D_NSN_LINESEARCH_enum")
      .value("SICONOS_FRICTION_3D_NSN_LINESEARCH_GOLDSTEINPRICE",
             SICONOS_FRICTION_3D_NSN_LINESEARCH_GOLDSTEINPRICE)
      .value("SICONOS_FRICTION_3D_NSN_LINESEARCH_ARMIJO",
             SICONOS_FRICTION_3D_NSN_LINESEARCH_ARMIJO)
      .value("SICONOS_FRICTION_3D_NSN_LINESEARCH_NO", SICONOS_FRICTION_3D_NSN_LINESEARCH_NO)
      .export_values();

  py::enum_<SICONOS_FRICTION_3D_NSN_HYBRID>(params, "SICONOS_FRICTION_3D_NSN_HYBRID_enum")
      .value("SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY_NO",
             SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY_NO)
      .value("SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY_PLI_NSN_LOOP",
             SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY_PLI_NSN_LOOP)
      .value("SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY_NSN_AND_PLI_NSN_LOOP",
             SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY_NSN_AND_PLI_NSN_LOOP)
      .value("SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY_VI_EG_NSN",
             SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY_VI_EG_NSN)
      .export_values();

  py::enum_<SICONOS_FRICTION_3D_PROXIMAL_IPARAM>(params,
                                                 "SICONOS_FRICTION_3D_PROXIMAL_IPARAM_enum")
      .value("SICONOS_FRICTION_3D_FP_ERROR_STRATEGY", SICONOS_FRICTION_3D_FP_ERROR_STRATEGY)
      .value("SICONOS_FRICTION_3D_PROXIMAL_IPARAM_CUMULATIVE_ITER_DONE",
             SICONOS_FRICTION_3D_PROXIMAL_IPARAM_CUMULATIVE_ITER_DONE)
      .value("SICONOS_FRICTION_3D_PROXIMAL_IPARAM_RELAXATION",
             SICONOS_FRICTION_3D_PROXIMAL_IPARAM_RELAXATION)
      .value("SICONOS_FRICTION_3D_PROXIMAL_IPARAM_STRATEGY",
             SICONOS_FRICTION_3D_PROXIMAL_IPARAM_STRATEGY)
      .export_values();

  py::enum_<SICONOS_FRICTION_3D_PROXIMAL_DPARAM>(params,
                                                 "SICONOS_FRICTION_3D_PROXIMAL_DPARAM_enum")
      .value("SICONOS_FRICTION_3D_PROXIMAL_DPARAM_ALPHA",
             SICONOS_FRICTION_3D_PROXIMAL_DPARAM_ALPHA)
      .value("SICONOS_FRICTION_3D_PROXIMAL_DPARAM_SIGMA",
             SICONOS_FRICTION_3D_PROXIMAL_DPARAM_SIGMA)
      .value("SICONOS_FRICTION_3D_PROXIMAL_DPARAM_NU", SICONOS_FRICTION_3D_PROXIMAL_DPARAM_NU)
      .value("SICONOS_FRICTION_3D_PROXIMAL_DPARAM_RELAXATION",
             SICONOS_FRICTION_3D_PROXIMAL_DPARAM_RELAXATION)
      .export_values();

  py::enum_<SICONOS_FRICTION_3D_PROXIMAL>(params, "SICONOS_FRICTION_3D_PROXIMAL_enum")
      .value("SICONOS_FRICTION_3D_PROXIMAL_PROX", SICONOS_FRICTION_3D_PROXIMAL_PROX)
      .value("SICONOS_FRICTION_3D_PROXIMAL_REGULARIZATION",
             SICONOS_FRICTION_3D_PROXIMAL_REGULARIZATION)
      .export_values();

  py::enum_<SICONOS_FRICTION_3D_ADMM_IPARAM>(params, "SICONOS_FRICTION_3D_ADMM_IPARAM_enum")
      .value("SICONOS_FRICTION_3D_ADMM_IPARAM_RHO_STRATEGY",
             SICONOS_FRICTION_3D_ADMM_IPARAM_RHO_STRATEGY)
      .value("SICONOS_FRICTION_3D_ADMM_IPARAM_INITIAL_RHO",
             SICONOS_FRICTION_3D_ADMM_IPARAM_INITIAL_RHO)
      .value("SICONOS_FRICTION_3D_ADMM_IPARAM_ACCELERATION",
             SICONOS_FRICTION_3D_ADMM_IPARAM_ACCELERATION)
      .value("SICONOS_FRICTION_3D_ADMM_IPARAM_SYMMETRY",
             SICONOS_FRICTION_3D_ADMM_IPARAM_SYMMETRY)
      .value("SICONOS_FRICTION_3D_ADMM_IPARAM_SPARSE_STORAGE",
             SICONOS_FRICTION_3D_ADMM_IPARAM_SPARSE_STORAGE)
      .value("SICONOS_FRICTION_3D_ADMM_IPARAM_GET_PROBLEM_INFO",
             SICONOS_FRICTION_3D_ADMM_IPARAM_GET_PROBLEM_INFO)
      .value("SICONOS_FRICTION_3D_ADMM_IPARAM_UPDATE_S",
             SICONOS_FRICTION_3D_ADMM_IPARAM_UPDATE_S)
      .value("SICONOS_FRICTION_3D_ADMM_IPARAM_FULL_H", SICONOS_FRICTION_3D_ADMM_IPARAM_FULL_H)
      .export_values();

  py::enum_<SICONOS_FRICTION_3D_ADMM_DPARAM>(params, "SICONOS_FRICTION_3D_ADMM_DPARAM_enum")
      .value("SICONOS_FRICTION_3D_ADMM_RHO", SICONOS_FRICTION_3D_ADMM_RHO)
      .value("SICONOS_FRICTION_3D_ADMM_RESTART_ETA", SICONOS_FRICTION_3D_ADMM_RESTART_ETA)
      .value("SICONOS_FRICTION_3D_ADMM_BALANCING_RESIDUAL_TAU",
             SICONOS_FRICTION_3D_ADMM_BALANCING_RESIDUAL_TAU)
      .value("SICONOS_FRICTION_3D_ADMM_BALANCING_RESIDUAL_PHI",
             SICONOS_FRICTION_3D_ADMM_BALANCING_RESIDUAL_PHI)
      .export_values();

  py::enum_<enum SICONOS_FRICTION_3D_ADMM_ACCELERATION>(
      params, "SICONOS_FRICTION_3D_ADMM_ACCELERATION")
      .value("SICONOS_FRICTION_3D_ADMM_NO_ACCELERATION_enum",
             SICONOS_FRICTION_3D_ADMM_NO_ACCELERATION)
      .value("SICONOS_FRICTION_3D_ADMM_ACCELERATION", SICONOS_FRICTION_3D_ADMM_ACCELERATION)
      .value("SICONOS_FRICTION_3D_ADMM_ACCELERATION_AND_RESTART",
             SICONOS_FRICTION_3D_ADMM_ACCELERATION_AND_RESTART)
      .export_values();

  py::enum_<SICONOS_FRICTION_3D_ADMM_SYMMETRY>(params,
                                               "SICONOS_FRICTION_3D_ADMM_SYMMETRY_enum")
      .value("SICONOS_FRICTION_3D_ADMM_CHECK_SYMMETRY",
             SICONOS_FRICTION_3D_ADMM_CHECK_SYMMETRY)
      .value("SICONOS_FRICTION_3D_ADMM_FORCED_SYMMETRY",
             SICONOS_FRICTION_3D_ADMM_FORCED_SYMMETRY)
      .value("SICONOS_FRICTION_3D_ADMM_FORCED_ASYMMETRY",
             SICONOS_FRICTION_3D_ADMM_FORCED_ASYMMETRY)
      .value("SICONOS_FRICTION_3D_ADMM_SYMMETRIZE", SICONOS_FRICTION_3D_ADMM_SYMMETRIZE)
      .value("SICONOS_FRICTION_3D_ADMM_ASSUME_SYMMETRY",
             SICONOS_FRICTION_3D_ADMM_ASSUME_SYMMETRY)
      .export_values();

  py::enum_<SICONOS_FRICTION_3D_ADMM_STORAGE>(params, "SICONOS_FRICTION_3D_ADMM_STORAGE_enum")
      .value("SICONOS_FRICTION_3D_ADMM_KEEP_STORAGE", SICONOS_FRICTION_3D_ADMM_KEEP_STORAGE)
      .value("SICONOS_FRICTION_3D_ADMM_FORCED_SPARSE_STORAGE",
             SICONOS_FRICTION_3D_ADMM_FORCED_SPARSE_STORAGE)
      .export_values();

  py::enum_<SICONOS_FRICTION_3D_ADMM_GET_PROBLEM_INFO>(
      params, "SICONOS_FRICTION_3D_ADMM_GET_PROBLEM_INFO_enum")
      .value("SICONOS_FRICTION_3D_ADMM_GET_PROBLEM_INFO_NO",
             SICONOS_FRICTION_3D_ADMM_GET_PROBLEM_INFO_NO)
      .value("SICONOS_FRICTION_3D_ADMM_GET_PROBLEM_INFO_YES",
             SICONOS_FRICTION_3D_ADMM_GET_PROBLEM_INFO_YES)
      .export_values();

  py::enum_<SICONOS_FRICTION_3D_ADMM_UPDATE_S>(params,
                                               "SICONOS_FRICTION_3D_ADMM_UPDATE_S_enum")
      .value("SICONOS_FRICTION_3D_ADMM_UPDATE_S_YES", SICONOS_FRICTION_3D_ADMM_UPDATE_S_YES)
      .value("SICONOS_FRICTION_3D_ADMM_UPDATE_S_NO", SICONOS_FRICTION_3D_ADMM_UPDATE_S_NO)
      .export_values();

  py::enum_<SICONOS_FRICTION_3D_ADMM_FULL_H>(params, "SICONOS_FRICTION_3D_ADMM_FULL_H_enum")
      .value("SICONOS_FRICTION_3D_ADMM_FULL_H_NO", SICONOS_FRICTION_3D_ADMM_FULL_H_NO)
      .value("SICONOS_FRICTION_3D_ADMM_FULL_H_YES", SICONOS_FRICTION_3D_ADMM_FULL_H_YES)
      .export_values();

  py::enum_<SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY>(
      params, "SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_enum")
      .value("SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_CONSTANT",
             SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_CONSTANT)
      .value("SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_RESIDUAL_BALANCING",
             SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_RESIDUAL_BALANCING)
      .value("SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_SCALED_RESIDUAL_BALANCING",
             SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_SCALED_RESIDUAL_BALANCING)
      .export_values();

  py::enum_<SICONOS_FRICTION_3D_ADMM_INITIAL_RHO>(params,
                                                  "SICONOS_FRICTION_3D_ADMM_INITIAL_RHO_enum")
      .value("SICONOS_FRICTION_3D_ADMM_INITIAL_RHO_GIVEN",
             SICONOS_FRICTION_3D_ADMM_INITIAL_RHO_GIVEN)
      .value("SICONOS_FRICTION_3D_ADMM_INITIAL_RHO_NORM_INF",
             SICONOS_FRICTION_3D_ADMM_INITIAL_RHO_NORM_INF)
      .value("SICONOS_FRICTION_3D_ADMM_INITIAL_RHO_EIGENVALUES",
             SICONOS_FRICTION_3D_ADMM_INITIAL_RHO_EIGENVALUES)
      .export_values();

  // Exposer l'énumération SICONOS_FRICTION_3D_IPM_IPARAM
  py::enum_<SICONOS_FRICTION_3D_IPM_IPARAM>(params, "SICONOS_FRICTION_3D_IPM_IPARAM")
      .value("SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_enum",
             SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING)
      .value("SICONOS_FRICTION_3D_IPM_IPARAM_SPARSE_STORAGE",
             SICONOS_FRICTION_3D_IPM_IPARAM_SPARSE_STORAGE)
      .value("SICONOS_FRICTION_3D_IPM_IPARAM_GET_PROBLEM_INFO",
             SICONOS_FRICTION_3D_IPM_IPARAM_GET_PROBLEM_INFO)
      .value("SICONOS_FRICTION_3D_IPM_IPARAM_ITERATES_MATLAB_FILE",
             SICONOS_FRICTION_3D_IPM_IPARAM_ITERATES_MATLAB_FILE)
      .value("SICONOS_FRICTION_3D_IPM_IPARAM_REDUCED_SYSTEM",
             SICONOS_FRICTION_3D_IPM_IPARAM_REDUCED_SYSTEM)
      .value("SICONOS_FRICTION_3D_IPM_IPARAM_FINISH_WITHOUT_SCALING",
             SICONOS_FRICTION_3D_IPM_IPARAM_FINISH_WITHOUT_SCALING)
      .value("SICONOS_FRICTION_3D_IPM_IPARAM_UPDATE_S",
             SICONOS_FRICTION_3D_IPM_IPARAM_UPDATE_S)
      .value("SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD",
             SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD)
      .value("SICONOS_FRICTION_3D_IPM_IPARAM_REDUCED_SYSTEM_METHOD",
             SICONOS_FRICTION_3D_IPM_IPARAM_REDUCED_SYSTEM_METHOD)
      .value("SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM", SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM)
      .value("SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT",
             SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT)
      .value("SICONOS_FRICTION_3D_IPM_IPARAM_CHOLESKY",
             SICONOS_FRICTION_3D_IPM_IPARAM_CHOLESKY)
      .value("SICONOS_FRICTION_3D_IPM_IPARAM_ITERATES_PYTHON_FILE",
             SICONOS_FRICTION_3D_IPM_IPARAM_ITERATES_PYTHON_FILE)
      .export_values();

  // Exposer l'énumération SICONOS_FRICTION_3D_IPM_DPARAM
  py::enum_<SICONOS_FRICTION_3D_IPM_DPARAM>(params, "SICONOS_FRICTION_3D_IPM_DPARAM_enum")
      .value("SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_1",
             SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_1)
      .value("SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_2",
             SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_2)
      .value("SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_3",
             SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_3)
      .value("SICONOS_FRICTION_3D_IPM_GAMMA_PARAMETER_1",
             SICONOS_FRICTION_3D_IPM_GAMMA_PARAMETER_1)
      .value("SICONOS_FRICTION_3D_IPM_GAMMA_PARAMETER_2",
             SICONOS_FRICTION_3D_IPM_GAMMA_PARAMETER_2)
      .export_values();

  // Exposer l'énumération SICONOS_FRICTION_3D_IPM_STORAGE
  py::enum_<SICONOS_FRICTION_3D_IPM_STORAGE>(params, "SICONOS_FRICTION_3D_IPM_STORAGE_enum")
      .value("SICONOS_FRICTION_3D_IPM_KEEP_STORAGE", SICONOS_FRICTION_3D_IPM_KEEP_STORAGE)
      .value("SICONOS_FRICTION_3D_IPM_FORCED_SPARSE_STORAGE",
             SICONOS_FRICTION_3D_IPM_FORCED_SPARSE_STORAGE)
      .export_values();

  // Exposer l'énumération SICONOS_FRICTION_3D_IPM_GET_PROBLEM_INFO
  py::enum_<SICONOS_FRICTION_3D_IPM_GET_PROBLEM_INFO>(
      params, "SICONOS_FRICTION_3D_IPM_GET_PROBLEM_INFO_enum")
      .value("SICONOS_FRICTION_3D_IPM_GET_PROBLEM_INFO_NO",
             SICONOS_FRICTION_3D_IPM_GET_PROBLEM_INFO_NO)
      .value("SICONOS_FRICTION_3D_IPM_GET_PROBLEM_INFO_YES",
             SICONOS_FRICTION_3D_IPM_GET_PROBLEM_INFO_YES)
      .export_values();

  // Exposer l'énumération SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_METHOD
  py::enum_<SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_METHOD>(
      params, "SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_METHOD_enum")
      .value("SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_QP",
             SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_QP)
      .value("SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_F",
             SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_F)
      .export_values();

  // Exposer l'énumération SICONOS_FRICTION_3D_IPM_IPARAM_REDUCED_SYSTEM_METHOD
  py::enum_<enum SICONOS_FRICTION_3D_IPM_IPARAM_REDUCED_SYSTEM_METHOD>(
      params, "SICONOS_FRICTION_3D_IPM_IPARAM_REDUCED_SYSTEM_METHOD_enum")
      .value("SICONOS_FRICTION_3D_IPM_IPARAM_REDUCED_SYSTEM_WITH_QP2",
             SICONOS_FRICTION_3D_IPM_IPARAM_REDUCED_SYSTEM_WITH_QP2)
      .value("SICONOS_FRICTION_3D_IPM_IPARAM_REDUCED_SYSTEM_WITH_QPH",
             SICONOS_FRICTION_3D_IPM_IPARAM_REDUCED_SYSTEM_WITH_QPH)
      .export_values();

  // Exposer l'énumération SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM
  py::enum_<enum SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM>(
      params, "SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM_enum")
      .value("SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_NOSCAL",
             SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_NOSCAL)
      .value("SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_QP2",
             SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_QP2)
      .value("SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_QPH",
             SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_QPH)
      .value("SICONOS_FRICTION_3D_IPM_IPARAM_LS_2X2_QP2",
             SICONOS_FRICTION_3D_IPM_IPARAM_LS_2X2_QP2)
      .value("SICONOS_FRICTION_3D_IPM_IPARAM_LS_2X2_QPH",
             SICONOS_FRICTION_3D_IPM_IPARAM_LS_2X2_QPH)
      .value("SICONOS_FRICTION_3D_IPM_IPARAM_LS_1X1_QPH",
             SICONOS_FRICTION_3D_IPM_IPARAM_LS_1X1_QPH)
      .value("SICONOS_FRICTION_3D_IPM_IPARAM_LS_1X1_JQJ",
             SICONOS_FRICTION_3D_IPM_IPARAM_LS_1X1_JQJ)
      .value("SICONOS_FRICTION_3D_IPM_IPARAM_LS_2X2_JQJ",
             SICONOS_FRICTION_3D_IPM_IPARAM_LS_2X2_JQJ)
      .value("SICONOS_FRICTION_3D_IPM_IPARAM_LS_2X2_invPH",
             SICONOS_FRICTION_3D_IPM_IPARAM_LS_2X2_invPH)
      .value("SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_JQinv",
             SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_JQinv)
      .export_values();

  // Exposer l'énumération SICONOS_FRICTION_3D_IPM_IPARAM_FINISH_WITHOUT_SCALING
  py::enum_<enum SICONOS_FRICTION_3D_IPM_IPARAM_FINISH_WITHOUT_SCALING>(
      params, "SICONOS_FRICTION_3D_IPM_IPARAM_FINISH_WITHOUT_SCALING_enum")
      .value("SICONOS_FRICTION_3D_IPM_IPARAM_FINISH_WITHOUT_SCALING_NO",
             SICONOS_FRICTION_3D_IPM_IPARAM_FINISH_WITHOUT_SCALING_NO)
      .value("SICONOS_FRICTION_3D_IPM_IPARAM_FINISH_WITHOUT_SCALING_YES",
             SICONOS_FRICTION_3D_IPM_IPARAM_FINISH_WITHOUT_SCALING_YES)
      .export_values();

  // Exposer l'énumération SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT
  py::enum_<enum SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT>(
      params, "SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT_enum")
      .value("SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT_NO",
             SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT_NO)
      .value("SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT_YES",
             SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT_YES)
      .value("SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT_AFTER",
             SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT_AFTER)
      .export_values();

  // Exposer l'énumération SICONOS_FRICTION_3D_IPM_IPARAM_CHOLESKY
  py::enum_<enum SICONOS_FRICTION_3D_IPM_IPARAM_CHOLESKY>(
      params, "SICONOS_FRICTION_3D_IPM_IPARAM_CHOLESKY_enum")
      .value("SICONOS_FRICTION_3D_IPM_IPARAM_CHOLESKY_NO",
             SICONOS_FRICTION_3D_IPM_IPARAM_CHOLESKY_NO)
      .value("SICONOS_FRICTION_3D_IPM_IPARAM_CHOLESKY_YES",
             SICONOS_FRICTION_3D_IPM_IPARAM_CHOLESKY_YES)
      .export_values();
}
