/* Siconos-Kernel, Copyright INRIA 2005-2015
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

/*! \file ControlTypeDef.hpp
  \brief Typedef for control-related objects
  */

#ifndef ControlTypeDef_H
#define ControlTypeDef_H

/** Actuator types */
#define PID_                       100
#define LINEAR_SMC                 101
#define EXPLICIT_LINEAR_SMC        103
#define LINEAR_SMC_OT2             104
#define LINEAR_SMC_IMPROVED        105
#define TWISTING                   106


/** Sensor types */
#define LINEAR_SENSOR              100

/** Observer types */
#define LUENBERGER                 100
#define SLIDING_REDUCED_ORDER      101

/** Event types
  \warning You have also to update the
  description in Event.hpp
*/
#define SENSOR_EVENT               3
#define OBSERVER_EVENT             4
#define ACTUATOR_EVENT             5

/** Base type forward declaration */

#endif
