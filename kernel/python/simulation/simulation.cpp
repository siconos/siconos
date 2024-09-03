#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <memory>

#include "TimeDiscretisation.hpp"
#include "TimeStepping.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "OneStepIntegrator.hpp"
#include "OneStepNSProblem.hpp"

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

PYBIND11_MODULE(simulation, m)
{
    m.doc() = "Siconos simulation module";

    py::class_<siconos::simulation::TimeDiscretisation, std::shared_ptr<siconos::simulation::TimeDiscretisation>>(m, "TimeDiscretisation")
        .def(py::init<double, double>());

    py::class_<siconos::simulation::Simulation, std::shared_ptr<siconos::simulation::Simulation>>(m, "Simulation")
        .def("hasNextEvent", &siconos::simulation::Simulation::hasNextEvent)
        .def("nextTime", &siconos::simulation::Simulation::nextTime);

    py::class_<siconos::simulation::TimeStepping, std::shared_ptr<siconos::simulation::TimeStepping>, siconos::simulation::Simulation>(m, "TimeStepping")
        .def(py::init<std::shared_ptr<siconos::modeling::NonSmoothDynamicalSystem>,
               std::shared_ptr<siconos::simulation::TimeDiscretisation>,
               std::shared_ptr<siconos::integrators::OneStepIntegrator>,
               std::shared_ptr<siconos::nonsmooth_formulations::OneStepNSProblem>>())
        .def("computeOneStep", &siconos::simulation::TimeStepping::computeOneStep)
        .def("nextStep", &siconos::simulation::TimeStepping::nextStep);
}