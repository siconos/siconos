#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "LagrangianLinearTIDS.hpp"
#include "SiconosVector.hpp"
#include "SimpleMatrix.hpp"

namespace py = pybind11;

void init_dynamical_systems(py::module &m)
{
  // LagrangianDS
  py::class_<siconos::modeling::LagrangianDS,
             std::shared_ptr<siconos::modeling::LagrangianDS>>(m, "LagrangianDS")
      .def(py::init<std::shared_ptr<siconos::algebra::SiconosVector>,
                    std::shared_ptr<siconos::algebra::SiconosVector>,
                    std::shared_ptr<siconos::algebra::SiconosMatrix>>(),
           py::arg().none(false), py::arg().none(false), py::arg().none(false))
      .def("forces", &siconos::modeling::LagrangianDS::forces)
      .def("computeForces", &siconos::modeling::LagrangianDS::computeForces)
      .def_property("fext", &siconos::modeling::LagrangianDS::fExt,
                    &siconos::modeling::LagrangianDS::setFExtPtr)
      .def("__repr__", [](const siconos::modeling::LagrangianDS &a) {
        a.display();
        return "\n";
      });

  // LagrangianLinearTIDS
  py::class_<siconos::modeling::LagrangianLinearTIDS, siconos::modeling::LagrangianDS,
             std::shared_ptr<siconos::modeling::LagrangianLinearTIDS>>(
      m, "LagrangianLinearTIDS", py::buffer_protocol())
      .def(py::init<std::shared_ptr<siconos::algebra::SiconosVector>,
                    std::shared_ptr<siconos::algebra::SiconosVector>,
                    std::shared_ptr<siconos::algebra::SiconosMatrix>>(),
           py::arg().none(false), py::arg().none(false), py::arg().none(false))
      .def("__repr__",
           [](const siconos::modeling::LagrangianLinearTIDS &a) {
             a.display();
             return "\n";
           })
      .def(py::init([](py::array q0, std::shared_ptr<siconos::algebra::SiconosVector> v0,
                       std::shared_ptr<siconos::algebra::SiconosMatrix> mass) {
        /* Request a buffer descriptor from Python */
        py::buffer_info info = q0.request();
        if (info.ndim != 1) throw std::runtime_error("Incompatible buffer dimension!");
        py::print("info ptr ", info.ptr);
        py::print("ref_count ", q0.ref_count());
        py::print(q0);
	//py::print(q0.data);
	std::vector<double> tmp = {};
	// tmp.data= std::move(q0.data);
	//py::print("tmp", tmp);
        // m.cast<siconos::algebra::SiconosVector>>(q0);
        //   auto qq0 = q0.cast<std::shared_ptr<siconos::algebra::SiconosVector>>();
        //   auto vv0 = v0.cast<std::shared_ptr<siconos::algebra::SiconosVector>>();
        //   auto mass0 = mass.cast<std::shared_ptr<siconos::algebra::SimpleMatrix>>();
        //   return siconos::modeling::LagrangianLinearTIDS(qq0, vv0, mass0);
        return std::make_shared<siconos::modeling::LagrangianLinearTIDS>(v0, v0, mass);
      }));

}
