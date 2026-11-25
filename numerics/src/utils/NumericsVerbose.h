/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */

/*!\file NumericsVerbose.h
 * \brief verbosity related functions (public) */

#ifndef _NUMERICS_VERBOSE_H_
#define _NUMERICS_VERBOSE_H_

/** logging types */
enum numerics_loggers {
  NUMERICS_LOG_TO_SCREEN,
  NUMERICS_LOG_TO_FILE,
  NUMERICS_EXTERNAL_LOGGER
};

#if defined(__cplusplus)
extern "C" {
#endif

/* Set verbose mode in numerics
   \param newVerboseMode 0 no verbose, 1 verbose.
*/
void numerics_set_verbose(int newVerboseMode);

#if defined(__cplusplus)
}
#endif

#endif
