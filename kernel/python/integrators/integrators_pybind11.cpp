#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <memory>

#include "MoreauJeanOSI.hpp"
#include "EulerMoreauOSI.hpp"

namespace py = pybind11;

// void init_integrators(py::module &m)
// {
//   // OSI base class
//   py::class_<siconos::integrators::OneStepIntegrator,
//              std::shared_ptr<siconos::integrators::OneStepIntegrator>>(m, "OneStepIntegrator");

//   // MoreauJean
//   py::class_<siconos::integrators::MoreauJeanOSI, siconos::integrators::OneStepIntegrator,
//              std::shared_ptr<siconos::integrators::MoreauJeanOSI>>(m, "MoreauJeanOSI")
//       // osi(theta, gamma)
//       .def(py::init<double, double>(), py::arg("theta") = 0.5,
//            py::arg("gamma") = std::numeric_limits<double>::quiet_NaN())

//       .def_property("theta", &siconos::integrators::MoreauJeanOSI::theta,
//                     &siconos::integrators::MoreauJeanOSI::setTheta)
//       .def_property("gamma", &siconos::integrators::MoreauJeanOSI::gamma,
//                     &siconos::integrators::MoreauJeanOSI::setGamma)

//       .def("__repr__", [](const siconos::integrators::MoreauJeanOSI &a) {
//         a.display();
//         return "\n";
//       });

//   // EulerMoreau
//   py::class_<siconos::integrators::EulerMoreauOSI, siconos::integrators::OneStepIntegrator,
//              std::shared_ptr<siconos::integrators::EulerMoreauOSI>>(m, "EulerMoreauOSI")
//       // osi(theta, gamma)
//       .def(py::init<double, double>(), py::arg("theta") = 0.5,
//            py::arg("gamma") = std::numeric_limits<double>::quiet_NaN())

//       .def_property("theta", &siconos::integrators::EulerMoreauOSI::theta,
//                     &siconos::integrators::EulerMoreauOSI::setTheta)
//       .def_property("gamma", &siconos::integrators::EulerMoreauOSI::gamma,
//                     &siconos::integrators::EulerMoreauOSI::setGamma)
//       .def_property("useGamma", &siconos::integrators::EulerMoreauOSI::useGamma,
//                     &siconos::integrators::EulerMoreauOSI::setUseGamma)
//       .def_property("useGammaForRelation",
//                     &siconos::integrators::EulerMoreauOSI::useGammaForRelation,
//                     &siconos::integrators::EulerMoreauOSI::setUseGammaForRelation)

//       .def("__repr__", [](const siconos::integrators::EulerMoreauOSI &a) {
//         a.display();
//         return "\n";
//       });
// }

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

PYBIND11_MODULE(integrators, m)
{
  // Optional docstring
  m.doc() = "Siconos one-step integrators";

  py::class_<siconos::integrators::OneStepIntegrator, std::shared_ptr<siconos::integrators::OneStepIntegrator>>(m, "OneStepIntegrator");

  py::class_<siconos::integrators::MoreauJeanOSI, std::shared_ptr<siconos::integrators::MoreauJeanOSI>, siconos::integrators::OneStepIntegrator>(m, "MoreauJeanOSI")
    .def(py::init<double, double>(), py::arg("theta") = 0.5,
          py::arg("gamma") = std::numeric_limits<double>::quiet_NaN())
    
    // .def("tonche_mass", &siconos::integrators::MoreauJeanOSI::tonch_mass) 
    
    .def("__repr__", [](const siconos::integrators::MoreauJeanOSI &a) {
      a.display();
      return "\n";
    });

  py::class_<siconos::integrators::EulerMoreauOSI, std::shared_ptr<siconos::integrators::EulerMoreauOSI>, siconos::integrators::OneStepIntegrator>(m, "EulerMoreauOSI")
    .def(py::init<double>());
}
