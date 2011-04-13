/* Siconos-Kernel, Copyright INRIA 2005-2010.
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

/*! \file Debug.hpp
  \brief Some debug facilities (cf wikipedia)
*/

#ifndef DEBUG_HPP
#define DEBUG_HPP

#ifdef DEBUG_MESSAGES
#define DEBUG_WHERESTR  "%s:%d:\n"
#define DEBUG_WHEREARG  __FILE__, __LINE__
#define DEBUG_INTERNAL_PRINTF(...)       fprintf(stderr, __VA_ARGS__)
#define DEBUG_PRINTF(_fmt, ...)  DEBUG_INTERNAL_PRINTF(DEBUG_WHERESTR _fmt, DEBUG_WHEREARG, __VA_ARGS__)
#define DEBUG_PRINT( M )  DEBUG_PRINTF("%s",M)
#else
#define DEBUG_PRINTF(_fmt, ...)
#define DEBUG_PRINT(M)
#endif

#endif
