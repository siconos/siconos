/* Siconos-Numerics, Copyright INRIA 2005-2014.
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

#ifndef VI_CST_H
#define VI_CST_H
/** \file VI_cst.h */
/** \enum VI_SOLVER VI_cst.h
 * Enum that allows one to encode the list of solvers in a proper to avoid mispelling
 * with char * variables
 */
enum VI_SOLVER
{
  SICONOS_VI_EG = 1000,
  SICONOS_VI_FPP = 1001,
  SICONOS_VI_HP = 1002,
  SICONOS_VI_BOX_QI = 1003,
  SICONOS_VI_BOX_AVI = 1004
};

extern char *  SICONOS_VI_EG_STR ;
extern char *  SICONOS_VI_FPP_STR ;
extern char *  SICONOS_VI_HP_STR ;
extern char *  SICONOS_VI_BOX_QI_STR ;
extern char *  SICONOS_VI_BOX_AVI_STR ;

#endif
