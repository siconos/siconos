/* Siconos-Kernel, Copyright INRIA 2005-2012.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */

/*! \file DynamicalSystemTypes.hpp
  \brief enum of the available types for dynamical systems,
  plugin names ...
*/

#ifndef DYNAMICALSYSTEMTYPES_HPP
#define DYNAMICALSYSTEMTYPES_HPP
/* now use : Type::<Siconos class>
 * Type::value(*pointer) see SiconosVisitor.hpp
 * transformation was done with
(tags-query-replace "FONLDS" "Type::FirstOrderNonLinearDS")
(tags-query-replace "FOLDS" "Type::FirstOrderLinearDS")
(tags-query-replace "FOLTIDS" "Type::FirstOrderLinearTIDS")
(tags-query-replace "LNLDS" "Type::LagrangianDS")
(tags-query-replace "LLTIDS" "Type::LagrangianLinearTIDS")
(tags-query-replace "NENLDS" "Type::NewtonEulerDS")
(tags-query-replace "ds->getType()" "Type::value(*ds)")
(tags-query-replace "DS::TYPES dsType" "Type::Siconos dsType")
and some manual interventions ...
*/

// Plugged objects used in Dynamical systems
#include "PluginTypes.hpp"


#endif
