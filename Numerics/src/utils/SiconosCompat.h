/* Siconos-Numerics, Copyright INRIA 2005-2011.
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
#ifndef _SICONOS_COMPAT
#define _SICONOS_COMPAT

#if defined(_MSC_VER)
#include <float.h>
#define isnan(x) _isnan(x)
#define isinf(x) (!_finite(x) && !_isnan(x))
inline double fmax(double x, double y)
{
  if (isnan(y)) return x;
  else if (isnan(x)) return y;
  else
  {
    if (x > y) return x;
    else return y;
  }
}
inline double fmin(double x, double y)
{
  if (isnan(y)) return x;
  else if (isnan(x)) return y;
  else
  {
    if (x > y) return y;
    else return x;
  }
}
#define INFINITY (DBL_MAX+DBL_MAX)
#define NAN (INFINITY-INFINITY)
#define _USE_MATH_DEFINES
#define copysign _copysign

#endif
#endif

