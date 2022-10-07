#include <pybind11/pybind11.h>

namespace py = pybind11;

void init_dynamical_systems(py::module &);
void init_nonsmoothlaws(py::module &);
void init_relations(py::module &);

PYBIND11_MODULE(modeling, m)
{
  // Optional docstring
  m.doc() = "Siconos modeling library";

  init_dynamical_systems(m);
  init_nonsmoothlaws(m);
  init_relations(m);
}
