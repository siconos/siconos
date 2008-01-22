/* Siconos-Numerics version 2.1.1, Copyright INRIA 2005-2007.
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
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
 */

#ifndef Numerics_Options_H
#define Numerics_Options_H

/*!\file Numerics_Options.h
  \brief General options for Numerics functions, structures and so on (mainly used to send information from Kernel to Numerics).
  \author Franck Perignon
*/

/** Structure used to set general options of Numerics functions, structures and so on.
    \param verbose mode (0: off, 1: on)
*/
typedef struct
{
  int verbose;
} Numerics_Options;


/* Verbose mode */
extern int Verb;

#ifdef __cplusplus
extern "C" {
#endif

  /* Set global option for numerics
     \param opt, a Numerics_Options structure
   */
  void setNumericsOptions(Numerics_Options* opt);

#ifdef __cplusplus
}
#endif


#endif
