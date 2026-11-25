// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

/*! \file MatrixBaseAddons.hpp

  Eigen MatrixBase class extensions

  See https://eigen.tuxfamily.org/dox/TopicCustomizing_Plugins.html

*/
// Setters

inline void setValue(Index row, Index col, Scalar v) { this->operator()(row, col) = v; }
