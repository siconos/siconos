/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2025 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

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
      "rewriteVectorFromAbsoluteToBodyFrame",
      [](const siconos::algebra::SiconosVector7 &q,
         Eigen::Ref<siconos::algebra::SiconosVector3> v) {
        siconos::geometry::rewriteVectorFromAbsoluteToBodyFrame(q, v);
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
      "rewriteVectorFromBodyToAbsoluteFrame",
      [](const siconos::algebra::SiconosVector7 &q,
         Eigen::Ref<siconos::algebra::SiconosVector3> v) {
        siconos::geometry::rewriteVectorFromBodyToAbsoluteFrame(q, v);
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
