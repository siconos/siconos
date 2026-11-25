/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */

/*! \file EigenInclude.hpp

  Set files where Eigen extensions are declared/defined

  Warning: Must be included before any Eigen's header file.

  See https://eigen.tuxfamily.org/dox/TopicCustomizing_Plugins.html
*/
#ifndef EIGENINCLUDE_HPP
#define EIGENINCLUDE_HPP

#include <iomanip>   // IWYU pragma: keep - required for setprecision
#include <iostream>  // IWYU pragma: keep - TODO : get rid of this! Note FP: maybe we should use a pre-compiled file process?
// #include <memory>  // For shared_ptr in MatrixAddons // TODO : get rid of this!

// filename of plugin for extending the MatrixBase class
#define EIGEN_MATRIXBASE_PLUGIN "MatrixBaseAddons.hpp"

// filename of plugin for extending the PlainObjectBase class.
#define EIGEN_PLAINOBJECTBASE_PLUGIN "PlainObjectAddons.hpp"

// filename of plugin for extending the MapBase class.
#define EIGEN_MAPBASE_PLUGIN "MapBaseAddons.hpp"

#endif
