/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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

/*!\file sn_logger.h
 * \brief Common data structures used by loggers in Siconos/Numerics
 *
 * \author Olivier Huber
 */

#ifndef SN_LOGGER_H
#define SN_LOGGER_H

/**\enum SN_loglevels sn_logger.h
 * loglevels for the loggers in numerics*/
typedef enum {
  SN_LOGLEVEL_NO,
  SN_LOGLEVEL_BASIC,
  SN_LOGLEVEL_LIGHT,
  SN_LOGLEVEL_VEC,
  SN_LOGLEVEL_MAT,
  SN_LOGLEVEL_ALL
} SN_loglevels;

#define SN_LOG_LIGHT(log_lvl, expr) if (log_lvl >= SN_LOGLEVEL_LIGHT) expr;
#define SN_LOG_SCALAR(log_lvl, expr) SN_LOG_LIGHT(log_lvl, expr)
#define SN_LOG_VEC(log_lvl, expr) if (log_lvl >= SN_LOGLEVEL_VEC) expr;
#define SN_LOG_MAT(log_lvl, expr) if (log_lvl >= SN_LOGLEVEL_MAT) expr;

#endif
