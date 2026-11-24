// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

/*!\file sn_logger.h
 * \brief Common data structures used by loggers in Siconos/Numerics
 *
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

#define SN_LOG_LIGHT(log_lvl, expr) \
  if (log_lvl >= SN_LOGLEVEL_LIGHT) expr;
#define SN_LOG_SCALAR(log_lvl, expr) SN_LOG_LIGHT(log_lvl, expr)
#define SN_LOG_VEC(log_lvl, expr) \
  if (log_lvl >= SN_LOGLEVEL_VEC) expr;
#define SN_LOG_MAT(log_lvl, expr) \
  if (log_lvl >= SN_LOGLEVEL_MAT) expr;

#endif
