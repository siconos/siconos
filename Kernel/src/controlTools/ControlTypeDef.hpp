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

/*! \file

Typedef for control-related objects
*/

#ifndef ControlTypeDef_H
#define ControlTypedef_H

/** Actuator types */
#define SAMPLED_PID_ACTUATOR  100
#define LINEAR_SMC            101
#define LINEAR_CHATTERING_SMC 103
#define LINEAR_SMC_OT2        104

/** Sensor types */
#define LINEAR_SENSOR         100

/** Event types
\warning Event.hpp has also to be updated
*/
#define SENSOR_EVENT          3
#define ACTUATOR_EVENT        4

#endif
