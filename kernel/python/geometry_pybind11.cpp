// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "RotationQuaternion.hpp"
#include "SiconosVector.hpp"

namespace py = pybind11;

PYBIND11_MODULE(geometry, m) {
  m.doc() = "Siconos tools to deal with quaternions, rotations ...";

  m.def("quaternionFromRotationVector", &siconos::geometry::quaternionFromRotationVector,
        py::arg("rotationVector"),
        R"pbdoc(
            Computes a quaternion from a rotation vector
            Parameters
            ----------
            rotationVector : array-like
               rotation (axis-angle)
            Returns
            -------
            numpy.ndarray of size 7
        )pbdoc");

  m.def("rotationVectorFromQuaternion", &siconos::geometry::rotationVectorFromQuaternion,
        py::arg("q0"), py::arg("q1"), py::arg("q2"), py::arg("q3"),
        R"pbdoc(
            Computes a quaternion from a rotation vector
            Parameters
            ----------
              q0, q1, q2, q3 quaternion components
            Returns
            -------
            numpy.ndarray of size 3
        )pbdoc");

  m.def(
      "rotateVectorFromInertialToBodyFrame",
      [](const siconos::algebra::SiconosVector7& q,
         Eigen::Ref<siconos::algebra::SiconosVector3> v) {
        siconos::geometry::rotateVectorFromInertialToBodyFrame(q, v);
      },
      py::arg("q"), py::arg("v"),
      R"pbdoc(
          For a given  configuration vector q composed of a position and a quaternion,
          express the vector v given in the inertial frame into to the body frame
          w.r.t the quaternion that parametrize the rotation in q.

          Parameters
          ----------
          q : array(7), [x,y,z,qw,qx,qy,qz]
          v : array(3)
        )pbdoc");

  m.def(
      "rotateVectorFromBodyToInertialFrame",
      [](const siconos::algebra::SiconosVector7& q,
         Eigen::Ref<siconos::algebra::SiconosVector3> v) {
        siconos::geometry::rotateVectorFromBodyToInertialFrame(q, v);
      },
      py::arg("q"), py::arg("v"),
      R"pbdoc(
          For a given  configuration vector q composed of a position and a quaternion,
          write an input vector expressed in the inertial frame into to the body frame
          w.r.t the quaternion that parametrize the rotation in q.
          Parameters
          ----------
          q : array(7), [x,y,z,qw,qx,qy,qz]
          v : array(3)
        )pbdoc");
}
