// # Siconos is a program dedicated to modeling, simulation and control
// # of non smooth dynamical systems.
// #
// # Copyright 2024 INRIA.
// #
// # Licensed under the Apache License, Version 2.0 (the "License");
// # you may not use this file except in compliance with the License.
// # You may obtain a copy of the License at
// #
// # http://www.apache.org/licenses/LICENSE-2.0
// #
// # Unless required by applicable law or agreed to in writing, software
// # distributed under the License is distributed on an "AS IS" BASIS,
// # WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// # See the License for the specific language governing permissions and
// # limitations under the License.
// #
// #


#include "Friction_cst.h"
#include "SolverOptions.h"
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>


namespace py = pybind11;

// Getter pour iparam qui retourne un tableau numpy lié directement au pointeur C++
py::array_t<int> get_iparam(SolverOptions &options) {
    return py::array_t<int>(
        {options.iSize},       // Dimension du tableau (taille de iparam)
        {sizeof(int)},         // Stride (décalage entre éléments)
        options.iparam,        // Pointeur vers les données C++ (iparam)
        py::cast(&options));   // Lier le cycle de vie de l'objet Python à l'objet C++
}


// Getter pour dparam qui retourne un tableau numpy lié directement au pointeur C++
py::array_t<double> get_dparam(SolverOptions &options) {
    return py::array_t<double>(
        {options.dSize},       // Dimension du tableau (taille de dparam)
        {sizeof(double)},      // Stride (décalage entre éléments)
        options.dparam,        // Pointeur vers les données C++ (dparam)
        py::cast(&options));   // Lier le cycle de vie de l'objet Python à l'objet C++
}


PYBIND11_MODULE(numerics, m) {

    py::class_<SolverOptions, std::shared_ptr<SolverOptions>>(m, "SolverOptions")
        // Membres simples
        .def(py::init<>())  // Constructeur par défaut
        .def_readwrite("solverId", &SolverOptions::solverId)
        .def_readwrite("isSet", &SolverOptions::isSet)
        .def_property_readonly("iparam", &get_iparam)  // Exposer iparam avec accès direct
        .def_property_readonly("dparam", &get_dparam)  // Exposer dparam avec accès direct
        .def_readwrite("filterOn", &SolverOptions::filterOn)
        .def("print", [](SolverOptions &options) {
            solver_options_print(&options);  // Utiliser la fonction existante pour afficher
        })
        .def("__repr__", [](const SolverOptions &options) {
            std::ostringstream oss;
            oss << "SolverOptions (ID: " << options.solverId << ")";
            solver_options_print(const_cast<SolverOptions*>(&options));
            return oss.str();
        });
        
        
    m.def("solver_options_create", &solver_options_create, py::return_value_policy::take_ownership, py::arg("solverId"));

    // Créer une classe "Constants" dans Python pour regrouper les enums
    py::class_<py::object>(m, "Constants")
        // SICONOS_IPARAM enum
        .def_property_readonly_static("SICONOS_IPARAM_MAX_ITER", [](py::object) { return static_cast<int>(SICONOS_IPARAM::SICONOS_IPARAM_MAX_ITER); })
        .def_property_readonly_static("SICONOS_IPARAM_ITER_DONE", [](py::object) { return static_cast<int>(SICONOS_IPARAM::SICONOS_IPARAM_ITER_DONE); })
        .def_property_readonly_static("SICONOS_IPARAM_PREALLOC", [](py::object) { return static_cast<int>(SICONOS_IPARAM::SICONOS_IPARAM_PREALLOC); })
        .def_property_readonly_static("SICONOS_IPARAM_NSGS_SHUFFLE", [](py::object) { return static_cast<int>(SICONOS_IPARAM::SICONOS_IPARAM_NSGS_SHUFFLE); })
        .def_property_readonly_static("SICONOS_IPARAM_ERROR_EVALUATION", [](py::object) { return static_cast<int>(SICONOS_IPARAM::SICONOS_IPARAM_ERROR_EVALUATION); })
        .def_property_readonly_static("SICONOS_IPARAM_PATHSEARCH_STACKSIZE", [](py::object) { return static_cast<int>(SICONOS_IPARAM::SICONOS_IPARAM_PATHSEARCH_STACKSIZE); })
        
        
        // SICONOS_DPARAM enum
        .def_property_readonly_static("SICONOS_DPARAM_TOL", [](py::object) { return static_cast<int>(SICONOS_DPARAM::SICONOS_DPARAM_TOL); })
        .def_property_readonly_static("SICONOS_DPARAM_RESIDU", [](py::object) { return static_cast<int>(SICONOS_DPARAM::SICONOS_DPARAM_RESIDU); })
        
        // FRICTION_SOLVER enum
        .def_property_readonly_static("SICONOS_FRICTION_2D_NSGS", [](py::object) { return static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_2D_NSGS); })
        .def_property_readonly_static("SICONOS_FRICTION_2D_CPG", [](py::object) { return static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_2D_CPG); })
        .def_property_readonly_static("SICONOS_FRICTION_2D_LEMKE", [](py::object) { return static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_2D_LEMKE); })
        .def_property_readonly_static("SICONOS_FRICTION_2D_ENUM", [](py::object) { return static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_2D_ENUM); })
        .def_property_readonly_static("SICONOS_FRICTION_3D_NSGS", [](py::object) { return static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_3D_NSGS); })
        .def_property_readonly_static("SICONOS_FRICTION_3D_NSGSV", [](py::object) { return static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_3D_NSGSV); })
        .def_property_readonly_static("SICONOS_FRICTION_3D_PROX", [](py::object) { return static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_3D_PROX); })
        .def_property_readonly_static("SICONOS_FRICTION_3D_TFP", [](py::object) { return static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_3D_TFP); })
        .def_property_readonly_static("SICONOS_FRICTION_3D_NSN_AC", [](py::object) { return static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_3D_NSN_AC); })
        .def_property_readonly_static("SICONOS_FRICTION_3D_DSFP", [](py::object) { return static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_3D_DSFP); })
        .def_property_readonly_static("SICONOS_FRICTION_3D_VI_FPP", [](py::object) { return static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_3D_VI_FPP); })
        .def_property_readonly_static("SICONOS_FRICTION_3D_VI_EG", [](py::object) { return static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_3D_VI_EG); })
        .def_property_readonly_static("SICONOS_FRICTION_3D_HP", [](py::object) { return static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_3D_HP); })
        .def_property_readonly_static("SICONOS_FRICTION_3D_NCPGlockerFBFixedPoint", [](py::object) { return static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_3D_NCPGlockerFBFixedPoint); })
        .def_property_readonly_static("SICONOS_FRICTION_3D_FPP", [](py::object) { return static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_3D_FPP); })
        .def_property_readonly_static("SICONOS_FRICTION_3D_EG", [](py::object) { return static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_3D_EG); })
        .def_property_readonly_static("SICONOS_FRICTION_3D_NSN_FB", [](py::object) { return static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_3D_NSN_FB); })
        .def_property_readonly_static("SICONOS_FRICTION_3D_GAMS_PATH", [](py::object) { return static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_3D_GAMS_PATH); })
        .def_property_readonly_static("SICONOS_FRICTION_3D_GAMS_PATHVI", [](py::object) { return static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_3D_GAMS_PATHVI); })
        .def_property_readonly_static("SICONOS_FRICTION_3D_ACLMFP", [](py::object) { return static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_3D_ACLMFP); })
        .def_property_readonly_static("SICONOS_FRICTION_3D_SOCLCP", [](py::object) { return static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_3D_SOCLCP); })
        .def_property_readonly_static("SICONOS_FRICTION_3D_GAMS_LCP_PATH", [](py::object) { return static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_3D_GAMS_LCP_PATH); })
        .def_property_readonly_static("SICONOS_FRICTION_3D_GAMS_LCP_PATHVI", [](py::object) { return static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_3D_GAMS_LCP_PATHVI); })
        .def_property_readonly_static("SICONOS_FRICTION_3D_NSN_NM", [](py::object) { return static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_3D_NSN_NM); })
        .def_property_readonly_static("SICONOS_FRICTION_3D_NSN_AC_TEST", [](py::object) { return static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_3D_NSN_AC_TEST); })
        .def_property_readonly_static("SICONOS_FRICTION_3D_PFP", [](py::object) { return static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_3D_PFP); })
        .def_property_readonly_static("SICONOS_FRICTION_3D_ADMM", [](py::object) { return static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_3D_ADMM); })
        .def_property_readonly_static("SICONOS_FRICTION_3D_ONECONTACT_NSN", [](py::object) { return static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_3D_ONECONTACT_NSN); })
        .def_property_readonly_static("SICONOS_FRICTION_3D_ONECONTACT_NSN_GP", [](py::object) { return static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_3D_ONECONTACT_NSN_GP); })
        .def_property_readonly_static("SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCone", [](py::object) { return static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCone); })
        .def_property_readonly_static("SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration", [](py::object) { return static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration); })
        .def_property_readonly_static("SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithRegularization", [](py::object) { return static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithRegularization); })
        .def_property_readonly_static("SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithDiagonalization", [](py::object) { return static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithDiagonalization); })
        .def_property_readonly_static("SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCone_velocity", [](py::object) { return static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCone_velocity); })
        .def_property_readonly_static("SICONOS_FRICTION_3D_NCPGlockerFBPATH", [](py::object) { return static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_3D_NCPGlockerFBPATH); })
        .def_property_readonly_static("SICONOS_FRICTION_3D_NCPGlockerFBNewton", [](py::object) { return static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_3D_NCPGlockerFBNewton); })
        .def_property_readonly_static("SICONOS_FRICTION_3D_ONECONTACT_QUARTIC", [](py::object) { return static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_3D_ONECONTACT_QUARTIC); })
        .def_property_readonly_static("SICONOS_FRICTION_3D_ONECONTACT_QUARTIC_NU", [](py::object) { return static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_3D_ONECONTACT_QUARTIC_NU); })
        .def_property_readonly_static("SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCylinder", [](py::object) { return static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCylinder); })
        .def_property_readonly_static("SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCylinderWithLocalIteration", [](py::object) { return static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCylinderWithLocalIteration); })
        .def_property_readonly_static("SICONOS_FRICTION_3D_ONECONTACT_NSN_GP_HYBRID", [](py::object) { return static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_3D_ONECONTACT_NSN_GP_HYBRID); })
        .def_property_readonly_static("SICONOS_FRICTION_3D_VI_FPP_Cylinder", [](py::object) { return static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_3D_VI_FPP_Cylinder); })
        .def_property_readonly_static("SICONOS_FRICTION_3D_CONVEXQP_PG_CYLINDER", [](py::object) { return static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_3D_CONVEXQP_PG_CYLINDER); })
        .def_property_readonly_static("SICONOS_GLOBAL_FRICTION_3D_NSGS_WR", [](py::object) { return static_cast<int>(FRICTION_SOLVER::SICONOS_GLOBAL_FRICTION_3D_NSGS_WR); })
        .def_property_readonly_static("SICONOS_GLOBAL_FRICTION_3D_NSGSV_WR", [](py::object) { return static_cast<int>(FRICTION_SOLVER::SICONOS_GLOBAL_FRICTION_3D_NSGSV_WR); })
        .def_property_readonly_static("SICONOS_GLOBAL_FRICTION_3D_PROX_WR", [](py::object) { return static_cast<int>(FRICTION_SOLVER::SICONOS_GLOBAL_FRICTION_3D_PROX_WR); })
        .def_property_readonly_static("SICONOS_GLOBAL_FRICTION_3D_DSFP_WR", [](py::object) { return static_cast<int>(FRICTION_SOLVER::SICONOS_GLOBAL_FRICTION_3D_DSFP_WR); })
        .def_property_readonly_static("SICONOS_GLOBAL_FRICTION_3D_TFP_WR", [](py::object) { return static_cast<int>(FRICTION_SOLVER::SICONOS_GLOBAL_FRICTION_3D_TFP_WR); })
        .def_property_readonly_static("SICONOS_GLOBAL_FRICTION_3D_NSGS", [](py::object) { return static_cast<int>(FRICTION_SOLVER::SICONOS_GLOBAL_FRICTION_3D_NSGS); })
        .def_property_readonly_static("SICONOS_GLOBAL_FRICTION_3D_NSN_AC_WR", [](py::object) { return static_cast<int>(FRICTION_SOLVER::SICONOS_GLOBAL_FRICTION_3D_NSN_AC_WR); })
        .def_property_readonly_static("SICONOS_GLOBAL_FRICTION_3D_NSN_AC", [](py::object) { return static_cast<int>(FRICTION_SOLVER::SICONOS_GLOBAL_FRICTION_3D_NSN_AC); })
        .def_property_readonly_static("SICONOS_GLOBAL_FRICTION_3D_GAMS_PATH", [](py::object) { return static_cast<int>(FRICTION_SOLVER::SICONOS_GLOBAL_FRICTION_3D_GAMS_PATH); })
        .def_property_readonly_static("SICONOS_GLOBAL_FRICTION_3D_GAMS_PATHVI", [](py::object) { return static_cast<int>(FRICTION_SOLVER::SICONOS_GLOBAL_FRICTION_3D_GAMS_PATHVI); })
        .def_property_readonly_static("SICONOS_GLOBAL_FRICTION_3D_VI_FPP", [](py::object) { return static_cast<int>(FRICTION_SOLVER::SICONOS_GLOBAL_FRICTION_3D_VI_FPP); })
        .def_property_readonly_static("SICONOS_GLOBAL_FRICTION_3D_VI_EG", [](py::object) { return static_cast<int>(FRICTION_SOLVER::SICONOS_GLOBAL_FRICTION_3D_VI_EG); })
        .def_property_readonly_static("SICONOS_GLOBAL_FRICTION_3D_ACLMFP", [](py::object) { return static_cast<int>(FRICTION_SOLVER::SICONOS_GLOBAL_FRICTION_3D_ACLMFP); })
        .def_property_readonly_static("SICONOS_GLOBAL_FRICTION_3D_ADMM", [](py::object) { return static_cast<int>(FRICTION_SOLVER::SICONOS_GLOBAL_FRICTION_3D_ADMM); })
        .def_property_readonly_static("SICONOS_GLOBAL_FRICTION_3D_ADMM_WR", [](py::object) { return static_cast<int>(FRICTION_SOLVER::SICONOS_GLOBAL_FRICTION_3D_ADMM_WR); })
        .def_property_readonly_static("SICONOS_GLOBAL_FRICTION_3D_IPM", [](py::object) { return static_cast<int>(FRICTION_SOLVER::SICONOS_GLOBAL_FRICTION_3D_IPM); })
        .def_property_readonly_static("SICONOS_ROLLING_FRICTION_3D_NSGS", [](py::object) { return static_cast<int>(FRICTION_SOLVER::SICONOS_ROLLING_FRICTION_3D_NSGS); })
        .def_property_readonly_static("SICONOS_ROLLING_FRICTION_3D_ONECONTACT_ProjectionOnCone", [](py::object) { return static_cast<int>(FRICTION_SOLVER::SICONOS_ROLLING_FRICTION_3D_ONECONTACT_ProjectionOnCone); })
        .def_property_readonly_static("SICONOS_ROLLING_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration", [](py::object) { return static_cast<int>(FRICTION_SOLVER::SICONOS_ROLLING_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration); })
        .def_property_readonly_static("SICONOS_ROLLING_FRICTION_3D_ADMM", [](py::object) { return static_cast<int>(FRICTION_SOLVER::SICONOS_ROLLING_FRICTION_3D_ADMM); })
        .def_property_readonly_static("SICONOS_ROLLING_FRICTION_2D_NSGS", [](py::object) { return static_cast<int>(FRICTION_SOLVER::SICONOS_ROLLING_FRICTION_2D_NSGS); })
        .def_property_readonly_static("SICONOS_ROLLING_FRICTION_2D_ONECONTACT_ProjectionOnCone", [](py::object) { return static_cast<int>(FRICTION_SOLVER::SICONOS_ROLLING_FRICTION_2D_ONECONTACT_ProjectionOnCone); })
        .def_property_readonly_static("SICONOS_ROLLING_FRICTION_2D_ONECONTACT_ProjectionOnConeWithLocalIteration", [](py::object) { return static_cast<int>(FRICTION_SOLVER::SICONOS_ROLLING_FRICTION_2D_ONECONTACT_ProjectionOnConeWithLocalIteration); })
        .def_property_readonly_static("SICONOS_GLOBAL_ROLLING_FRICTION_3D_NSGS_WR", [](py::object) { return static_cast<int>(FRICTION_SOLVER::SICONOS_GLOBAL_ROLLING_FRICTION_3D_NSGS_WR); })
        .def_property_readonly_static("SICONOS_GLOBAL_ROLLING_FRICTION_3D_IPM", [](py::object) { return static_cast<int>(FRICTION_SOLVER::SICONOS_GLOBAL_ROLLING_FRICTION_3D_IPM); });

}