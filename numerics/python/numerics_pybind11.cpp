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
#include "relay_cst.h"
#include "lcp_cst.h"
#include "NM_types.h"
#include "GenericMechanical_cst.h"
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


PYBIND11_MODULE(pynumerics, m) {

    py::class_<SolverOptions, std::shared_ptr<SolverOptions>>(m, "SolverOptions")
        // Membres simples
        .def(py::init<>())  // Constructeur par défaut
        .def_readwrite("solverId", &SolverOptions::solverId)
        .def_readonly("iSize", &SolverOptions::iSize)
        .def_readonly("dSize", &SolverOptions::dSize)
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

    py::module_ constants = m.def_submodule("constants", "Constants for numerics module");

    py::enum_<NumericsMatrix_types>(constants, "NumericsMatrix_types", "Types of storage for NumericsMatrix")
    .value("NM_DENSE", NM_DENSE, "Dense format")
    .value("NM_SPARSE_BLOCK", NM_SPARSE_BLOCK, "Sparse block format")
    .value("NM_SPARSE", NM_SPARSE, "Compressed column format")
    .value("NM_UNKNOWN", NM_UNKNOWN, "Unset. Used in NM_null")
    .export_values();

    // SICONOS_IPARAM enum
    constants.attr("SICONOS_IPARAM_MAX_ITER") = static_cast<int>(SICONOS_IPARAM::SICONOS_IPARAM_MAX_ITER);
    constants.attr("SICONOS_IPARAM_ITER_DONE") = static_cast<int>(SICONOS_IPARAM::SICONOS_IPARAM_ITER_DONE);
    constants.attr("SICONOS_IPARAM_PREALLOC") = static_cast<int>(SICONOS_IPARAM::SICONOS_IPARAM_PREALLOC);
    constants.attr("SICONOS_IPARAM_NSGS_SHUFFLE") = static_cast<int>(SICONOS_IPARAM::SICONOS_IPARAM_NSGS_SHUFFLE);
    constants.attr("SICONOS_IPARAM_ERROR_EVALUATION") = static_cast<int>(SICONOS_IPARAM::SICONOS_IPARAM_ERROR_EVALUATION);
    constants.attr("SICONOS_IPARAM_PATHSEARCH_STACKSIZE") = static_cast<int>(SICONOS_IPARAM::SICONOS_IPARAM_PATHSEARCH_STACKSIZE);

    // SICONOS_DPARAM enum
    constants.attr("SICONOS_DPARAM_TOL") = static_cast<int>(SICONOS_DPARAM::SICONOS_DPARAM_TOL);
    constants.attr("SICONOS_DPARAM_RESIDU") = static_cast<int>(SICONOS_DPARAM::SICONOS_DPARAM_RESIDU);

    // FRICTION_SOLVER enum
    constants.attr("SICONOS_FRICTION_2D_NSGS") = static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_2D_NSGS);
    constants.attr("SICONOS_FRICTION_2D_CPG") = static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_2D_CPG);
    constants.attr("SICONOS_FRICTION_2D_LEMKE") = static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_2D_LEMKE);
    constants.attr("SICONOS_FRICTION_2D_ENUM") = static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_2D_ENUM);
    constants.attr("SICONOS_FRICTION_3D_NSGS") = static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_3D_NSGS);
    constants.attr("SICONOS_FRICTION_3D_NSGSV") = static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_3D_NSGSV);
    constants.attr("SICONOS_FRICTION_3D_PROX") = static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_3D_PROX);
    constants.attr("SICONOS_FRICTION_3D_TFP") = static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_3D_TFP);
    constants.attr("SICONOS_FRICTION_3D_NSN_AC") = static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_3D_NSN_AC);
    constants.attr("SICONOS_FRICTION_3D_DSFP") = static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_3D_DSFP);
    constants.attr("SICONOS_FRICTION_3D_VI_FPP") = static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_3D_VI_FPP);
    constants.attr("SICONOS_FRICTION_3D_VI_EG") = static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_3D_VI_EG);
    constants.attr("SICONOS_FRICTION_3D_HP") = static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_3D_HP);
    constants.attr("SICONOS_FRICTION_3D_NCPGlockerFBFixedPoint") = static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_3D_NCPGlockerFBFixedPoint);
    constants.attr("SICONOS_FRICTION_3D_FPP") = static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_3D_FPP);
    constants.attr("SICONOS_FRICTION_3D_EG") = static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_3D_EG);
    constants.attr("SICONOS_FRICTION_3D_NSN_FB") = static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_3D_NSN_FB);
    constants.attr("SICONOS_FRICTION_3D_GAMS_PATH") = static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_3D_GAMS_PATH);
    constants.attr("SICONOS_FRICTION_3D_GAMS_PATHVI") = static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_3D_GAMS_PATHVI);
    constants.attr("SICONOS_FRICTION_3D_ACLMFP") = static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_3D_ACLMFP);
    constants.attr("SICONOS_FRICTION_3D_SOCLCP") = static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_3D_SOCLCP);
    constants.attr("SICONOS_FRICTION_3D_GAMS_LCP_PATH") = static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_3D_GAMS_LCP_PATH);
    constants.attr("SICONOS_FRICTION_3D_GAMS_LCP_PATHVI") = static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_3D_GAMS_LCP_PATHVI);
    constants.attr("SICONOS_FRICTION_3D_NSN_NM") = static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_3D_NSN_NM);
    constants.attr("SICONOS_FRICTION_3D_NSN_AC_TEST") = static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_3D_NSN_AC_TEST);
    constants.attr("SICONOS_FRICTION_3D_PFP") = static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_3D_PFP);
    constants.attr("SICONOS_FRICTION_3D_ADMM") = static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_3D_ADMM);
    constants.attr("SICONOS_FRICTION_3D_ONECONTACT_NSN") = static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_3D_ONECONTACT_NSN);
    constants.attr("SICONOS_FRICTION_3D_ONECONTACT_NSN_GP") = static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_3D_ONECONTACT_NSN_GP);
    constants.attr("SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCone") = static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCone);
    constants.attr("SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration") = static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration);
    constants.attr("SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithRegularization") = static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithRegularization);
    constants.attr("SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithDiagonalization") = static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithDiagonalization);
    constants.attr("SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCone_velocity") = static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCone_velocity);
    constants.attr("SICONOS_FRICTION_3D_NCPGlockerFBPATH") = static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_3D_NCPGlockerFBPATH);
    constants.attr("SICONOS_FRICTION_3D_NCPGlockerFBNewton") = static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_3D_NCPGlockerFBNewton);
    constants.attr("SICONOS_FRICTION_3D_ONECONTACT_QUARTIC") = static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_3D_ONECONTACT_QUARTIC);
    constants.attr("SICONOS_FRICTION_3D_ONECONTACT_QUARTIC_NU") = static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_3D_ONECONTACT_QUARTIC_NU);
    constants.attr("SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCylinder") = static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCylinder);
    constants.attr("SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCylinderWithLocalIteration") = static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCylinderWithLocalIteration);
    constants.attr("SICONOS_FRICTION_3D_ONECONTACT_NSN_GP_HYBRID") = static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_3D_ONECONTACT_NSN_GP_HYBRID);
    constants.attr("SICONOS_FRICTION_3D_VI_FPP_Cylinder") = static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_3D_VI_FPP_Cylinder);
    constants.attr("SICONOS_FRICTION_3D_CONVEXQP_PG_CYLINDER") = static_cast<int>(FRICTION_SOLVER::SICONOS_FRICTION_3D_CONVEXQP_PG_CYLINDER);
    constants.attr("SICONOS_GLOBAL_FRICTION_3D_NSGS_WR") = static_cast<int>(FRICTION_SOLVER::SICONOS_GLOBAL_FRICTION_3D_NSGS_WR);
    constants.attr("SICONOS_GLOBAL_FRICTION_3D_NSGSV_WR") = static_cast<int>(FRICTION_SOLVER::SICONOS_GLOBAL_FRICTION_3D_NSGSV_WR);
    constants.attr("SICONOS_GLOBAL_FRICTION_3D_PROX_WR") = static_cast<int>(FRICTION_SOLVER::SICONOS_GLOBAL_FRICTION_3D_PROX_WR);
    constants.attr("SICONOS_GLOBAL_FRICTION_3D_DSFP_WR") = static_cast<int>(FRICTION_SOLVER::SICONOS_GLOBAL_FRICTION_3D_DSFP_WR);
    constants.attr("SICONOS_GLOBAL_FRICTION_3D_TFP_WR") = static_cast<int>(FRICTION_SOLVER::SICONOS_GLOBAL_FRICTION_3D_TFP_WR);
    constants.attr("SICONOS_GLOBAL_FRICTION_3D_NSGS") = static_cast<int>(FRICTION_SOLVER::SICONOS_GLOBAL_FRICTION_3D_NSGS);
    constants.attr("SICONOS_GLOBAL_FRICTION_3D_NSN_AC_WR") = static_cast<int>(FRICTION_SOLVER::SICONOS_GLOBAL_FRICTION_3D_NSN_AC_WR);
    constants.attr("SICONOS_GLOBAL_FRICTION_3D_NSN_AC") = static_cast<int>(FRICTION_SOLVER::SICONOS_GLOBAL_FRICTION_3D_NSN_AC);
    constants.attr("SICONOS_GLOBAL_FRICTION_3D_GAMS_PATH") = static_cast<int>(FRICTION_SOLVER::SICONOS_GLOBAL_FRICTION_3D_GAMS_PATH);
    constants.attr("SICONOS_GLOBAL_FRICTION_3D_GAMS_PATHVI") = static_cast<int>(FRICTION_SOLVER::SICONOS_GLOBAL_FRICTION_3D_GAMS_PATHVI);
    constants.attr("SICONOS_GLOBAL_FRICTION_3D_VI_FPP") = static_cast<int>(FRICTION_SOLVER::SICONOS_GLOBAL_FRICTION_3D_VI_FPP);
    constants.attr("SICONOS_GLOBAL_FRICTION_3D_VI_EG") = static_cast<int>(FRICTION_SOLVER::SICONOS_GLOBAL_FRICTION_3D_VI_EG);
    constants.attr("SICONOS_GLOBAL_FRICTION_3D_ACLMFP") = static_cast<int>(FRICTION_SOLVER::SICONOS_GLOBAL_FRICTION_3D_ACLMFP);
    constants.attr("SICONOS_GLOBAL_FRICTION_3D_ADMM") = static_cast<int>(FRICTION_SOLVER::SICONOS_GLOBAL_FRICTION_3D_ADMM);
    constants.attr("SICONOS_GLOBAL_FRICTION_3D_ADMM_WR") = static_cast<int>(FRICTION_SOLVER::SICONOS_GLOBAL_FRICTION_3D_ADMM_WR);
    constants.attr("SICONOS_GLOBAL_FRICTION_3D_IPM") = static_cast<int>(FRICTION_SOLVER::SICONOS_GLOBAL_FRICTION_3D_IPM);
    constants.attr("SICONOS_ROLLING_FRICTION_3D_NSGS") = static_cast<int>(FRICTION_SOLVER::SICONOS_ROLLING_FRICTION_3D_NSGS);
    constants.attr("SICONOS_ROLLING_FRICTION_3D_ONECONTACT_ProjectionOnCone") = static_cast<int>(FRICTION_SOLVER::SICONOS_ROLLING_FRICTION_3D_ONECONTACT_ProjectionOnCone);
    constants.attr("SICONOS_ROLLING_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration") = static_cast<int>(FRICTION_SOLVER::SICONOS_ROLLING_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration);
    constants.attr("SICONOS_ROLLING_FRICTION_3D_ADMM") = static_cast<int>(FRICTION_SOLVER::SICONOS_ROLLING_FRICTION_3D_ADMM);
    constants.attr("SICONOS_ROLLING_FRICTION_2D_NSGS") = static_cast<int>(FRICTION_SOLVER::SICONOS_ROLLING_FRICTION_2D_NSGS);
    constants.attr("SICONOS_ROLLING_FRICTION_2D_ONECONTACT_ProjectionOnCone") = static_cast<int>(FRICTION_SOLVER::SICONOS_ROLLING_FRICTION_2D_ONECONTACT_ProjectionOnCone);
    constants.attr("SICONOS_ROLLING_FRICTION_2D_ONECONTACT_ProjectionOnConeWithLocalIteration") = static_cast<int>(FRICTION_SOLVER::SICONOS_ROLLING_FRICTION_2D_ONECONTACT_ProjectionOnConeWithLocalIteration);
    constants.attr("SICONOS_GLOBAL_ROLLING_FRICTION_3D_NSGS_WR") = static_cast<int>(FRICTION_SOLVER::SICONOS_GLOBAL_ROLLING_FRICTION_3D_NSGS_WR);
    constants.attr("SICONOS_GLOBAL_ROLLING_FRICTION_3D_IPM") = static_cast<int>(FRICTION_SOLVER::SICONOS_GLOBAL_ROLLING_FRICTION_3D_IPM);

    // RELAY_SOLVER enum
    constants.attr("SICONOS_RELAY_PGS") = static_cast<int>(RELAY_SOLVER::SICONOS_RELAY_PGS);
    constants.attr("SICONOS_RELAY_ENUM") = static_cast<int>(RELAY_SOLVER::SICONOS_RELAY_ENUM);
    constants.attr("SICONOS_RELAY_PATH") = static_cast<int>(RELAY_SOLVER::SICONOS_RELAY_PATH);
    constants.attr("SICONOS_RELAY_LEMKE") = static_cast<int>(RELAY_SOLVER::SICONOS_RELAY_LEMKE);
    constants.attr("SICONOS_RELAY_AVI_CAOFERRIS") = static_cast<int>(RELAY_SOLVER::SICONOS_RELAY_AVI_CAOFERRIS);
    constants.attr("SICONOS_RELAY_AVI_CAOFERRIS_TEST") = static_cast<int>(RELAY_SOLVER::SICONOS_RELAY_AVI_CAOFERRIS_TEST);
    
    // LCP_SOLVER enum
    constants.attr("SICONOS_LCP_PGS") = static_cast<int>(LCP_SOLVER::SICONOS_LCP_PGS);
    constants.attr("SICONOS_LCP_ENUM") = static_cast<int>(LCP_SOLVER::SICONOS_LCP_ENUM);
    constants.attr("SICONOS_LCP_PATH") = static_cast<int>(LCP_SOLVER::SICONOS_LCP_PATH);
    constants.attr("SICONOS_LCP_LEMKE") = static_cast<int>(LCP_SOLVER::SICONOS_LCP_LEMKE);
    constants.attr("SICONOS_LCP_AVI_CAOFERRIS") = static_cast<int>(LCP_SOLVER::SICONOS_LCP_AVI_CAOFERRIS);
      
    // GENERIC_MECHANICAL_SOLVER enum
    constants.attr("SICONOS_GENERIC_MECHANICAL_NSGS") = static_cast<int>(GENERIC_MECHANICAL_SOLVER::SICONOS_GENERIC_MECHANICAL_NSGS);
}
