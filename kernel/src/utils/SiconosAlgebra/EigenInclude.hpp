/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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

/*! \file EigenInclude.hpp

  Set files where Eigen extensions are declared/defined

  Warning: Must be included before any Eigen's header file.

  See https://eigen.tuxfamily.org/dox/TopicCustomizing_Plugins.html
*/
#ifndef EIGENINCLUDE_HPP
#define EIGENINCLUDE_HPP

#include <iomanip>  // required for setprecision
#include <iostream>  // TODO : get rid of this! Note FP: maybe we should use a pre-compiled file process?
#include <memory>  // For shared_ptr in MatrixAddons // TODO : get rid of this!

// filename of plugin for extending the MatrixBase class
#define EIGEN_MATRIXBASE_PLUGIN "MatrixBaseAddons.hpp"

// filename of plugin for extending the PlainObjectBase class.
#define EIGEN_PLAINOBJECTBASE_PLUGIN "PlainObjectAddons.hpp"

// filename of plugin for extending the Matrix class.
#define EIGEN_MATRIX_PLUGIN "MatrixAddons.hpp"

// filename of plugin for extending the MapBase class.
#define EIGEN_MAPBASE_PLUGIN "MapBaseAddons.hpp"

#endif
