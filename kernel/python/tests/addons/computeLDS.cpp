/* cppimport

<%

import os
setup_pybind11(cfg)

with open("cppimport_header.py") as f:
    exec(f.read())

%>
*/
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <SiconosKernel.hpp>
#include <SiconosMatrix.hpp>
#include <SiconosVector.hpp>
#include <StorageTools.hpp>
#include <iostream>

namespace py = pybind11;

siconos::algebra::SiconosSparseMatrix initializeSparse(
    const Eigen::Ref<const siconos::algebra::SiconosVector>& q,
    siconos::algebra::SiconosSparseMatrix& mass) {
  mass.setZero();

  for (int i = 0; i < q.size(); ++i) mass.insert(i, i) = 1.0 + q(i) * 0.1;
  mass.makeCompressed();
  return mass;
}

siconos::algebra::SiconosSparseMatrix computeMassSparse(
    const Eigen::Ref<const siconos::algebra::SiconosVector>& q,
    siconos::algebra::SiconosSparseMatrix& mass) {
  for (int k = 0; k < mass.outerSize(); ++k) {
    for (siconos::algebra::SiconosSparseMatrix::InnerIterator it(mass, k); it; ++it) {
      it.valueRef() *= 1000.;
    }
  }
  return mass;
}

PYBIND11_MODULE(computeLDS, addons) {
  addons.def("initializeSparse", &initializeSparse);
  addons.def("computeMassSparse", &computeMassSparse);
}
