#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <memory>

#include "InteractionManager.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "OneStepIntegrator.hpp"
#include "OneStepNSProblem.hpp"
#include "SimulationGraphs.hpp"
#include "SimulationTypes.hpp"
#include "TimeDiscretisation.hpp"
#include "TimeStepping.hpp"
#include "TimeSteppingDirectProjection.hpp"
#include "Topology.hpp"

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

PYBIND11_MODULE(simulation, m) {
  m.doc() = "Siconos simulation module";

  py::class_<siconos::simulation::InteractionManager,
             std::shared_ptr<siconos::simulation::InteractionManager>>(m, "InteractionManager")
      .def("updateInteractions", &siconos::simulation::InteractionManager::updateInteractions)
      .def("insertNonSmoothLaw", &siconos::simulation::InteractionManager::insertNonSmoothLaw)
      .def("nonSmoothLaw", &siconos::simulation::InteractionManager::nonSmoothLaw);

  py::class_<siconos::simulation::TimeDiscretisation,
             std::shared_ptr<siconos::simulation::TimeDiscretisation>>(m, "TimeDiscretisation")
      .def(py::init<double, double>());

  py::class_<siconos::simulation::Topology, std::shared_ptr<siconos::simulation::Topology>>(
      m, "Topology")
      .def("indexSetsSize", &siconos::simulation::Topology::indexSetsSize)
      .def("numberOfIndexSet", &siconos::simulation::Topology::numberOfIndexSet)
      .def("getDynamicalSystem",
           py::overload_cast<std::string>(&siconos::simulation::Topology::getDynamicalSystem,
                                          py::const_))
      .def("getInteraction", py::overload_cast<std::string>(
                                 &siconos::simulation::Topology::getInteraction, py::const_));

  py::class_<siconos::simulation::Simulation,
             std::shared_ptr<siconos::simulation::Simulation>>(m, "Simulation")
      .def("hasNextEvent", &siconos::simulation::Simulation::hasNextEvent)
      .def("nextTime", &siconos::simulation::Simulation::nextTime)
      .def("insertNonSmoothProblem", &siconos::simulation::Simulation::insertNonSmoothProblem,
           py::arg("osns"), py::arg("Id") = 0, "Insert a non-smooth problem into the system")
      .def("insertInteractionManager",
           &siconos::simulation::Simulation::insertInteractionManager, py::arg("manager"),
           "Set an object to automatically manage interactions during the simulation")
      .def("startingTime", &siconos::simulation::Simulation::startingTime,
           "get current time (ie starting point for current integration, time of currentEvent "
           "of eventsManager.)")
      .def("oneStepNSProblem", &siconos::simulation::Simulation::oneStepNSProblem,
           "get the OneStepNSProblem with the given id")
      .def("clearNSDSChangeLog", &siconos::simulation::Simulation::clearNSDSChangeLog,
           "clear the change log of the non-smooth dynamical system")
      .def("initializeOSIAssociations",
           &siconos::simulation::Simulation::initializeOSIAssociations,
           "initialize integrators-DS links in the NSDS graph")
      .def("initializeIndexSets", &siconos::simulation::Simulation::initializeIndexSets,
           "initialize integrators index sets")
      .def("applyNSDSChangelogForDS",
           &siconos::simulation::Simulation::applyNSDSChangelogForDS,
           "initialize/update all new dynamical systems and interactions")
      .def("computeInitialStateOfTheStep",
           &siconos::simulation::Simulation::computeInitialStateOfTheStep)
      .def("updateWorldFromDS", &siconos::simulation::Simulation::updateWorldFromDS)
      .def("updateInteractions", &siconos::simulation::Simulation::updateInteractions,
           "Call the interaction manager")
      .def("initializeNSDSChangelog",
           &siconos::simulation::Simulation::initializeNSDSChangelog,
           " initialize objects (DSs and Interations) found in the NSDS")
      .def("updateOutput", &siconos::simulation::Simulation::updateOutput,
           py::arg("level") = 0)
      .def("updateIndexSets", &siconos::simulation::Simulation::updateIndexSets,
           " update all index sets of the topology, using current states of the"
           "Interactions")
      .def("initializeOneStepNSProblem",
           &siconos::simulation::Simulation::initializeOneStepNSProblem)
      .def("firstInitialize", &siconos::simulation::Simulation::firstInitialize)
      .def("updateDSPlugins", &siconos::simulation::Simulation::updateDSPlugins)
      .def("computeResidu", &siconos::simulation::Simulation::computeResidu,
           "Compute the residu for all integrators")
      .def("updateAllInput", &siconos::simulation::Simulation::updateAllInput)
      .def("updateState", &siconos::simulation::Simulation::updateState, py::arg("level") = 0)
      .def("y_output", &siconos::simulation::Simulation::y_output,
           py::return_value_policy::move)
      .def("lambda_input", &siconos::simulation::Simulation::lambda_input,
           py::return_value_policy::move)
      .def("indexSet", &siconos::simulation::Simulation::indexSet)
      .def_property_readonly("numberOfOSNSProblems",
                             &siconos::simulation::Simulation::numberOfOSNSProblems);

  py::class_<siconos::simulation::TimeStepping,
             std::shared_ptr<siconos::simulation::TimeStepping>,
             siconos::simulation::Simulation>(m, "TimeStepping")
      .def(py::init<std::shared_ptr<siconos::modeling::NonSmoothDynamicalSystem>,
                    std::shared_ptr<siconos::simulation::TimeDiscretisation>,
                    std::shared_ptr<siconos::integrators::OneStepIntegrator>,
                    std::shared_ptr<siconos::nonsmooth_formulations::OneStepNSProblem>>())
      .def(py::init<std::shared_ptr<siconos::modeling::NonSmoothDynamicalSystem>,
                    std::shared_ptr<siconos::simulation::TimeDiscretisation>, int>(),
           py::arg("nsds"), py::arg("td"), py::arg("nb") = 0)
      .def("computeOneStep", &siconos::simulation::TimeStepping::computeOneStep)
      .def("nextStep", &siconos::simulation::TimeStepping::nextStep)
      .def("timeStep", &siconos::simulation::TimeStepping::timeStep)
      .def("insertIntegrator", &siconos::simulation::TimeStepping::insertIntegrator)
      .def("setNewtonOptions", &siconos::simulation::TimeStepping::setNewtonOptions,
           py::arg("newton_options"))
      .def("setNewtonMaxIteration", &siconos::simulation::TimeStepping::setNewtonMaxIteration,
           py::arg("newton_max_iter"))
      .def("setNewtonTolerance", &siconos::simulation::TimeStepping::setNewtonTolerance,
           py::arg("newton_tol"))
      .def("setNewtonWarningOnNonConvergence",
           &siconos::simulation::TimeStepping::setNewtonWarningOnNonConvergence,
           py::arg("newton_warning"))
      .def("setWarningNonsmoothSolver",
           &siconos::simulation::TimeStepping::setWarningNonsmoothSolver,
           py::arg("nonsmooth_warning"))
      .def("setSkipLastUpdateOutput",
           &siconos::simulation::TimeStepping::setSkipLastUpdateOutput,
           py::arg("skip_last_update_output"))
      .def("setSkipLastUpdateInput",
           &siconos::simulation::TimeStepping::setSkipLastUpdateInput,
           py::arg("skip_last_update_input"))
      .def("setSkipResetLambdas", &siconos::simulation::TimeStepping::setSkipResetLambdas,
           py::arg("skip_reset_lambdas"))
      .def_property_readonly("skipResetLambdas",
                             &siconos::simulation::TimeStepping::skipResetLambdas)
      .def("setDisplayNewtonConvergence",
           &siconos::simulation::TimeStepping::setDisplayNewtonConvergence,
           py::arg("display_newton_convergence"))
      .def("resetLambdas", &siconos::simulation::TimeStepping::resetLambdas)
      .def_property_readonly("newtonTolerance",
                             &siconos::simulation::TimeStepping::newtonTolerance)
      .def_property_readonly("newtonMaxIteration",
                             &siconos::simulation::TimeStepping::newtonMaxIteration)
      .def_property_readonly("newtonOptions",
                             &siconos::simulation::TimeStepping::newtonOptions)
      .def("computeFreeState", &siconos::simulation::TimeStepping::computeFreeState)
      .def("DefaultCheckSolverOutput",
           &siconos::simulation::TimeStepping::DefaultCheckSolverOutput)
      .def("skipLastUpdateInput", &siconos::simulation::TimeStepping::skipLastUpdateInput)
      .def("skipLastUpdateOutput", &siconos::simulation::TimeStepping::skipLastUpdateOutput)
      .def("prepareNewtonIteration",
           &siconos::simulation::TimeStepping::prepareNewtonIteration);

  py::class_<siconos::simulation::TimeSteppingDirectProjection,
             std::shared_ptr<siconos::simulation::TimeSteppingDirectProjection>,
             siconos::simulation::TimeStepping>(m, "TimeSteppingDirectProjection")
      .def(py::init<std::shared_ptr<siconos::modeling::NonSmoothDynamicalSystem>,
                    std::shared_ptr<siconos::simulation::TimeDiscretisation>,
                    std::shared_ptr<siconos::integrators::OneStepIntegrator>,
                    std::shared_ptr<siconos::nonsmooth_formulations::OneStepNSProblem>,
                    std::shared_ptr<siconos::nonsmooth_formulations::OneStepNSProblem>,
                    unsigned int>(),
           py::arg("nsds"), py::arg("td"), py::arg("osi"), py::arg("osnspb_velo"),
           py::arg("osnspb_pos"), py::arg("level") = 1)
      .def("advanceToEvent",
           &siconos::simulation::TimeSteppingDirectProjection::advanceToEvent)
      .def("nextStep", &siconos::simulation::TimeSteppingDirectProjection::nextStep)
      .def("computeCriteria",
           &siconos::simulation::TimeSteppingDirectProjection::computeCriteria)
      .def("setProjectionMaxIteration",
           &siconos::simulation::TimeSteppingDirectProjection::setProjectionMaxIteration,
           py::arg("v"))
      .def("setConstraintTolUnilateral",
           &siconos::simulation::TimeSteppingDirectProjection::setConstraintTolUnilateral,
           py::arg("v"))
      .def("setConstraintTol",
           &siconos::simulation::TimeSteppingDirectProjection::setConstraintTol, py::arg("v"));

  py::class_<siconos::graphs::InteractionsGraph,
             std::shared_ptr<siconos::graphs::InteractionsGraph>>(m, "InteractionsGraph");

  py::module_ constants = m.def_submodule("constants", "Constants for simulation module");
  constants.attr("SICONOS_OSNSP_TS_VELOCITY") =
      static_cast<int>(siconos::simulation::SICONOS_OSNSP_TS_VELOCITY);
  constants.attr("SICONOS_OSNSP_TS_POS") =
      static_cast<int>(siconos::simulation::SICONOS_OSNSP_TS_POS);

  py::enum_<siconos::simulation::TimeSteppingType>(m, "TimeSteppingType",
                                                   "define the type of time stepping scheme")
      .value("LINEAR", siconos::simulation::TimeSteppingType::LINEAR,
             "force a single iteration of the Newton Solver")
      .value("LINEAR_IMPLICIT", siconos::simulation::TimeSteppingType::LINEAR_IMPLICIT,
             "force a single iteration of the Newton Solver")
      .value("NONLINEAR", siconos::simulation::TimeSteppingType::NONLINEAR,
             "performs the newton iterations up to convergence")
      .value("NONLINEAR_FULL", siconos::simulation::TimeSteppingType::NONLINEAR_FULL,
             "performs the newton iterations up to convergence")
      .export_values();
}
