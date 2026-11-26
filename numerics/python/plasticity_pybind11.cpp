// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>

#include <cstring>
#include <memory>
#include <stdexcept>
#include <string>

#include "NonSmoothDrivers.h"
#include "NumericsMatrix.h"
#include "PlasticityProblem.h"
#include "Plasticity_options.h"
#include "SolverOptions.h"
#include "plasticity_2d_compute_error.h"

namespace py = pybind11;

namespace {
struct PlasticityProblemDeleter {
  void operator()(PlasticityProblem* problem) const {
    if (problem) plasticityProblem_free(problem);
  }
};
using PlasticityProblemPtr = std::unique_ptr<PlasticityProblem, PlasticityProblemDeleter>;

/* malloc'ed copy of a vector of the expected size */
double* copy_vector(const py::array_t<double>& v, ssize_t size, const char* name) {
  auto a = py::array_t<double, py::array::c_style | py::array::forcecast>::ensure(v);
  if (!a || a.size() != size)
    throw std::invalid_argument(std::string(name) + " must have size " + std::to_string(size));
  double* data = (double*)malloc((size_t)size * sizeof(double));
  std::memcpy(data, a.data(), (size_t)size * sizeof(double));
  return data;
}

/* problem with dimension, number of cones, M (dense, column-major copy) and q */
PlasticityProblemPtr new_problem(int dim, int nc, const py::array_t<double>& M,
                                 const py::array_t<double>& q) {
  ssize_t n = (ssize_t)dim * nc;
  auto Mf = py::array_t<double, py::array::f_style | py::array::forcecast>::ensure(M);
  if (!Mf || Mf.ndim() != 2 || Mf.shape(0) != n || Mf.shape(1) != n)
    throw std::invalid_argument("M must be a " + std::to_string(n) + "x" + std::to_string(n) +
                                " matrix");
  PlasticityProblemPtr problem(plasticityProblem_new());
  problem->dimension = dim;
  problem->numberOfCones = nc;
  double* data = (double*)malloc((size_t)(n * n) * sizeof(double));
  std::memcpy(data, Mf.data(), (size_t)(n * n) * sizeof(double));
  problem->M = NM_create_with_data(NM_DENSE, (int)n, (int)n, data);
  problem->q = copy_vector(q, n, "q");
  return problem;
}
}  // namespace

void wrap_plasticity(py::module_& m, py::module_& params, py::module_& solver_ids) {
  // PlasticityModelType enum
  py::enum_<PlasticityModelType>(m, "PlasticityModelType", "Type of plasticity model")
      .value("UNKNOWN", PLASTICITY_MODEL_UNKNOWN, "Unknown model")
      .value("DRUCKER_PRAGER", PLASTICITY_MODEL_DRUCKER_PRAGER, "Drucker-Prager model")
      .value("VON_MISES", PLASTICITY_MODEL_VON_MISES, "Von Mises model")
      .export_values();

  // PlasticityVonMisesHypothesis enum
  py::enum_<PlasticityVonMisesHypothesis>(m, "PlasticityVonMisesHypothesis",
                                          "Von Mises yield set")
      .value("PLANE_STRAIN", PLASTICITY_VONMISES_PLANE_STRAIN,
             "Plane strain, szz = (sxx + syy)/2: cylinder ||sigma_T|| <= sqrt(2/3) sigma_y")
      .value("PLANE_STRESS", PLASTICITY_VONMISES_PLANE_STRESS,
             "Plane stress, szz = 0: ellipsoid sigma_N^2/2 + 3/2 ||sigma_T||^2 <= sigma_y^2")
      .value("BALL", PLASTICITY_VONMISES_BALL, "Ball ||sigma|| <= sigma_y")
      .export_values();

  // PlasticityStressCoordinates enum
  py::enum_<PlasticityStressCoordinates>(m, "PlasticityStressCoordinates",
                                         "Coordinates of the stress (and strain rate)")
      .value("CONE", PLASTICITY_STRESS_COORDINATES_CONE,
             "Cone coordinates (sigma_N, sigma_T1, sigma_T2)")
      .value("CARTESIAN_MANDEL", PLASTICITY_STRESS_COORDINATES_CARTESIAN_MANDEL,
             "Cartesian, Mandel notation (sxx, syy, sqrt(2) sxy)")
      .value("CARTESIAN_VOIGT", PLASTICITY_STRESS_COORDINATES_CARTESIAN_VOIGT,
             "Cartesian, Voigt notation (sxx, syy, sxy), strain rate (exx, eyy, 2 exy)")
      .export_values();

  // Plasticity_DruckerPrager_model
  py::class_<Plasticity_DruckerPrager_model, py::smart_holder>(
      m, "Plasticity_DruckerPrager_model",
      "Drucker-Prager plasticity model parameters (eta, theta)")
      .def(py::init<>())
      .def_readwrite("eta", &Plasticity_DruckerPrager_model::eta, "Cone coefficients")
      .def_readwrite("theta", &Plasticity_DruckerPrager_model::theta, "Dilatancy coefficients")
      .def_readwrite("cohesion", &Plasticity_DruckerPrager_model::cohesion,
                     "Cohesions (may be None)");

  // Plasticity_VonMises_model
  py::class_<Plasticity_VonMises_model, py::smart_holder>(
      m, "Plasticity_VonMises_model", "Von Mises plasticity model parameters (yield stress)")
      .def(py::init<>())
      .def_readwrite("sigma_y", &Plasticity_VonMises_model::sigma_y, "Yield stress")
      .def_readwrite("hypothesis", &Plasticity_VonMises_model::hypothesis,
                     "Yield set: plane strain cylinder, plane stress ellipsoid or ball");

  // PlasticityProblem - use a wrapper approach similar to friction contact
  py::class_<PlasticityProblem, py::smart_holder>(m, "PlasticityProblem",
                                                  "Generic plasticity problem structure")
      .def(py::init([]() { return PlasticityProblemPtr(plasticityProblem_new()); }))
      .def_readwrite("dimension", &PlasticityProblem::dimension, "Dimension of stress space")
      .def_readwrite("numberOfCones", &PlasticityProblem::numberOfCones, "Number of cones")
      .def_readwrite("M", &PlasticityProblem::M, "M matrix")
      .def_readwrite("q", &PlasticityProblem::q, "q vector")
      .def_readwrite("model_type", &PlasticityProblem::model_type, "Model type enum")
      .def_readwrite("model", &PlasticityProblem::model, "Model parameters union")
      .def_readwrite("stress_coordinates", &PlasticityProblem::stress_coordinates,
                     "Coordinates of the stress and strain rate (cone or Cartesian)");

  // Factory functions. The problems are allocated by plasticityProblem_new and freed by
  // plasticityProblem_free (unique_ptr deleter kept by the smart holder). M is copied in
  // column-major order (NM_DENSE), whatever the order of the numpy array.
  m.def(
      "PlasticityProblem_new",
      [](int dim, py::array_t<double> M, py::array_t<double> q, py::array_t<double> eta,
         py::array_t<double> theta, py::object cohesion) {
        int nc = (int)eta.size();
        PlasticityProblemPtr problem = new_problem(dim, nc, M, q);
        problem->model_type = PLASTICITY_MODEL_DRUCKER_PRAGER;
        problem->model.drucker_prager =
            (Plasticity_DruckerPrager_model*)calloc(1, sizeof(Plasticity_DruckerPrager_model));
        problem->model.drucker_prager->eta = copy_vector(eta, nc, "eta");
        problem->model.drucker_prager->theta = copy_vector(theta, nc, "theta");
        // Optional cohesion (K = {|sigma_T| <= eta sigma_N + c})
        if (!cohesion.is_none())
          problem->model.drucker_prager->cohesion =
              copy_vector(cohesion.cast<py::array_t<double>>(), nc, "cohesion");
        return problem;
      },
      py::arg("dim"), py::arg("M"), py::arg("q"), py::arg("eta"), py::arg("theta"),
      py::arg("cohesion") = py::none(),
      "Create a PlasticityProblem with Drucker-Prager model, K = {|sigma_T| <= eta sigma_N + "
      "c}, cohesion c optional (handled by plasticity_2d_driver)");

  // Backward compatibility alias
  m.attr("Plasticity2DProblem") = m.attr("PlasticityProblem");

  // Backward compatibility factory function - same as PlasticityProblem_new
  m.def(
      "Plasticity2DProblem_new",
      [](int dim, py::array_t<double> M, py::array_t<double> q, py::array_t<double> eta,
         py::array_t<double> theta) {
        int nc = (int)eta.size();
        PlasticityProblemPtr problem = new_problem(dim, nc, M, q);
        problem->model_type = PLASTICITY_MODEL_DRUCKER_PRAGER;
        problem->model.drucker_prager =
            (Plasticity_DruckerPrager_model*)calloc(1, sizeof(Plasticity_DruckerPrager_model));
        problem->model.drucker_prager->eta = copy_vector(eta, nc, "eta");
        problem->model.drucker_prager->theta = copy_vector(theta, nc, "theta");
        return problem;
      },
      py::arg("dim"), py::arg("M"), py::arg("q"), py::arg("eta"), py::arg("theta"),
      "Create a PlasticityProblem (backward compatibility)");

  // Factory function for creating PlasticityProblem with Von Mises model
  m.def(
      "PlasticityProblem_new_VonMises",
      [](int dim, py::array_t<double> M, py::array_t<double> q, py::array_t<double> sigma_y,
         PlasticityVonMisesHypothesis hypothesis) {
        int nc = (int)sigma_y.size();
        PlasticityProblemPtr problem = new_problem(dim, nc, M, q);
        problem->model_type = PLASTICITY_MODEL_VON_MISES;
        problem->model.von_mises =
            (Plasticity_VonMises_model*)calloc(1, sizeof(Plasticity_VonMises_model));
        problem->model.von_mises->sigma_y = copy_vector(sigma_y, nc, "sigma_y");
        problem->model.von_mises->hypothesis = hypothesis;
        return problem;
      },
      py::arg("dim"), py::arg("M"), py::arg("q"), py::arg("sigma_y"),
      py::arg("hypothesis") = PLASTICITY_VONMISES_PLANE_STRAIN,
      "Create a PlasticityProblem with Von Mises model");

  // Driver: stress (initial guess / solution) and strainrate are modified in place
  m.def(
      "plasticity_2d_driver",
      [](PlasticityProblem* problem, py::array_t<double, py::array::c_style> stress,
         py::array_t<double, py::array::c_style> strainrate, SolverOptions* options) {
        ssize_t n = 3 * problem->numberOfCones;
        if (stress.size() != n || strainrate.size() != n)
          throw std::invalid_argument("stress and strainrate must have size 3 numberOfCones");
        SolverOptions* opts = options ? options : solver_options_create(PLASTICITY_2D_NSGS);
        int info = plasticity_2d_driver(problem, stress.mutable_data(),
                                        strainrate.mutable_data(), opts);
        if (!options) solver_options_delete(opts);
        return info;
      },
      py::arg("problem"), py::arg("stress"), py::arg("strainrate"),
      py::arg("options") = nullptr,
      "Solve the problem (default solver: PLASTICITY_2D_NSGS). stress (initial guess) and "
      "strainrate (float64, contiguous, size 3 numberOfCones) are modified in place. "
      "Returns info (0: converged).");

  // Error: strainrate is set to M stress + q
  m.def(
      "plasticity_2d_compute_error",
      [](PlasticityProblem* problem, py::array_t<double, py::array::c_style> stress,
         py::array_t<double, py::array::c_style> strainrate, double tolerance,
         SolverOptions* options) {
        ssize_t n = 3 * problem->numberOfCones;
        if (stress.size() != n || strainrate.size() != n)
          throw std::invalid_argument("stress and strainrate must have size 3 numberOfCones");
        double error = 0.0;
        int info = plasticity_2d_compute_error(problem, stress.mutable_data(),
                                               strainrate.mutable_data(), tolerance, options,
                                               0.0, &error);
        return py::make_tuple(info, error);
      },
      py::arg("problem"), py::arg("stress"), py::arg("strainrate"), py::arg("tolerance"),
      py::arg("options") = nullptr,
      "Natural-map error of the problem (cohesion, coordinates and Von Mises hypothesis taken "
      "into account), scaled according to options.iparam[params.PLASTICITY_SCALING_*] "
      "(default scaling if options is None); strainrate is set to M stress + q. Returns "
      "(info, error).");

  // Display functions
  m.def("plasticity_display", &plasticity_display, "Display plasticity problem",
        py::arg("problem"));

  // Backward compatibility function names
  m.def(
      "mohrCoulomb2D_display", [](PlasticityProblem* problem) { plasticity_display(problem); },
      "Display plasticity problem (backward compatibility)", py::arg("problem"));

  // Scaling of the error of plasticity_2d_compute_error (indices in iparam / dparam of the
  // NSGS options, same as the fc3d ones)
  params.attr("PLASTICITY_SCALING_NORM") = (int)PLASTICITY_SCALING_NORM;
  params.attr("PLASTICITY_SCALING_SCOPE") = (int)PLASTICITY_SCALING_SCOPE;
  params.attr("PLASTICITY_SCALING_REDUCTION") = (int)PLASTICITY_SCALING_REDUCTION;
  params.attr("PLASTICITY_SCALING_ALPHA") = (int)PLASTICITY_SCALING_ALPHA;
  params.attr("PLASTICITY_SCALING_RHO") = (int)PLASTICITY_SCALING_RHO;
  py::enum_<PLASTICITY_SCALING_RHO_ENUM>(params, "PLASTICITY_SCALING_RHO_enum",
                                         "rho of the natural map of the error")
      .value("PLASTICITY_SCALING_RHO_ONE", PLASTICITY_SCALING_RHO_ONE)
      .value("PLASTICITY_SCALING_RHO_ELEMENT", PLASTICITY_SCALING_RHO_ELEMENT)
      .value("PLASTICITY_SCALING_RHO_GLOBAL", PLASTICITY_SCALING_RHO_GLOBAL)
      .export_values();
  py::enum_<PLASTICITY_SCALING>(params, "PLASTICITY_SCALING", "Norm of the error scaling")
      .value("PLASTICITY_SCALING_EUCLIDEAN", PLASTICITY_SCALING_EUCLIDEAN)
      .value("PLASTICITY_SCALING_NO", PLASTICITY_SCALING_NO)
      .value("PLASTICITY_SCALING_CONE", PLASTICITY_SCALING_CONE)
      .value("PLASTICITY_SCALING_CONE_ETA", PLASTICITY_SCALING_CONE_ETA)
      .value("PLASTICITY_SCALING_CONE_WEIGHTED", PLASTICITY_SCALING_CONE_WEIGHTED)
      .export_values();
  py::enum_<enum PLASTICITY_SCALING_SCOPE>(params, "PLASTICITY_SCALING_SCOPE_enum",
                                           "Scope of the error scaling")
      .value("PLASTICITY_SCALING_SCOPE_PER_CONE", PLASTICITY_SCALING_SCOPE_PER_CONE)
      .value("PLASTICITY_SCALING_SCOPE_PER_VECTOR", PLASTICITY_SCALING_SCOPE_PER_VECTOR)
      .export_values();
  py::enum_<enum PLASTICITY_SCALING_REDUCTION>(params, "PLASTICITY_SCALING_REDUCTION_enum",
                                               "Reduction of the error scaling")
      .value("PLASTICITY_SCALING_REDUCTION_MAX", PLASTICITY_SCALING_REDUCTION_MAX)
      .value("PLASTICITY_SCALING_REDUCTION_MEAN", PLASTICITY_SCALING_REDUCTION_MEAN)
      .value("PLASTICITY_SCALING_REDUCTION_MIXED", PLASTICITY_SCALING_REDUCTION_MIXED)
      .export_values();

  // Solver IDs enum
  py::enum_<PLASTICITY_SOLVER>(solver_ids, "PLASTICITY_SOLVER", "Plasticity solver IDs")
      .value("PLASTICITY_2D_NSGS", PLASTICITY_2D_NSGS, "Non-smooth Gauss-Seidel")
      .value("PLASTICITY_2D_NSGS_GENERIC", PLASTICITY_2D_NSGS_GENERIC,
             "Non-smooth Gauss-Seidel (generic)")
      .value("PLASTICITY_2D_NSN", PLASTICITY_2D_NSN, "Global semi-smooth Newton")
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
      .value("PLASTICITY_2D_ONECONE_VONMISES_RETURN_MAPPING",
             PLASTICITY_2D_ONECONE_VONMISES_RETURN_MAPPING, "Von Mises return mapping (exact)")
      .value("PLASTICITY_2D_ONECONE_VONMISES_FP", PLASTICITY_2D_ONECONE_VONMISES_FP,
             "Von Mises fixed point on the natural map")
      .value("PLASTICITY_2D_ONECONE_VONMISES_NSN", PLASTICITY_2D_ONECONE_VONMISES_NSN,
             "Von Mises semi-smooth Newton on the natural map")
      .value("PLASTICITY_2D_ONECONE_DRUCKERPRAGER_RETURN_MAPPING",
             PLASTICITY_2D_ONECONE_DRUCKERPRAGER_RETURN_MAPPING,
             "Drucker-Prager, one cone: return mapping (exact enumeration)")
      .value("PLASTICITY_2D_ONECONE_DRUCKERPRAGER_NSN",
             PLASTICITY_2D_ONECONE_DRUCKERPRAGER_NSN,
             "Drucker-Prager, one cone: semi-smooth Newton on the natural map")
      .export_values();

  // Also export as module attributes for convenience
  m.attr("PLASTICITY_2D_NSGS") = (int)PLASTICITY_2D_NSGS;
  m.attr("PLASTICITY_2D_NSGS_GENERIC") = (int)PLASTICITY_2D_NSGS_GENERIC;
  m.attr("PLASTICITY_2D_NSN") = (int)PLASTICITY_2D_NSN;
  m.attr("PLASTICITY_2D_ONECONE_NSN") = (int)PLASTICITY_2D_ONECONE_NSN;
  m.attr("PLASTICITY_2D_ONECONE_NSN_GP") = (int)PLASTICITY_2D_ONECONE_NSN_GP;
  m.attr("PLASTICITY_2D_ONECONE_NSN_GP_HYBRID") = (int)PLASTICITY_2D_ONECONE_NSN_GP_HYBRID;
  m.attr("PLASTICITY_2D_ONECONE_ProjectionOnCone") =
      (int)PLASTICITY_2D_ONECONE_ProjectionOnCone;
  m.attr("PLASTICITY_2D_ONECONE_ProjectionOnConeWithLocalIteration") =
      (int)PLASTICITY_2D_ONECONE_ProjectionOnConeWithLocalIteration;
  m.attr("PLASTICITY_2D_ONECONE_VONMISES_RETURN_MAPPING") =
      (int)PLASTICITY_2D_ONECONE_VONMISES_RETURN_MAPPING;
  m.attr("PLASTICITY_2D_ONECONE_VONMISES_FP") = (int)PLASTICITY_2D_ONECONE_VONMISES_FP;
  m.attr("PLASTICITY_2D_ONECONE_VONMISES_NSN") = (int)PLASTICITY_2D_ONECONE_VONMISES_NSN;
  m.attr("PLASTICITY_2D_ONECONE_DRUCKERPRAGER_RETURN_MAPPING") =
      (int)PLASTICITY_2D_ONECONE_DRUCKERPRAGER_RETURN_MAPPING;
  m.attr("PLASTICITY_2D_ONECONE_DRUCKERPRAGER_NSN") =
      (int)PLASTICITY_2D_ONECONE_DRUCKERPRAGER_NSN;
  // deprecated name of PLASTICITY_2D_ONECONE_VONMISES_RETURN_MAPPING
  m.attr("PLASTICITY_2D_ONECONE_VONMISES") = (int)PLASTICITY_2D_ONECONE_VONMISES;

  // Backward compatibility solver IDs
  m.attr("MOHR_COULOMB_2D_NSGS") = (int)PLASTICITY_2D_NSGS;
  m.attr("SICONOS_MC2D_NSGS") = (int)PLASTICITY_2D_NSGS;
}
