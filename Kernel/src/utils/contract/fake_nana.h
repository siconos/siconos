/* Siconos-Kernel version 1.1.2, Copyright INRIA 2005-2006.
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
#ifndef FAKE_NANA_H
#define FAKE_NANA_H

/*
 * fake_nana.h -
 * provide header file when compile without nana
 * -> all the functions implemented in nana are redefined and do nothing
 */

#define VL(a) /* empty */
#define VLG(a) /* empty */
#define VLH(a) /* empty */
#define VLP(a) /* empty */
#define VLGP(a) /* empty */
#define VLHP(a) /* empty */
#define VLGHP(a) /* empty */
#define I(e) /* empty */
#define IG(e,g) /* empty */
#define IH(e,h) /* empty */
#define IP(e,p) /* empty */
#define IGH(e,g,h) /* empty */
#define IGP(e,g,p) /* empty */
#define IHP(e,h,p) /* empty */
#define IGHP(e,g,h,p) /* empty */

#define N(e) /* empty */
#define NG(e,g) /* empty */
#define NH(e,h) /* empty */
#define NP(e,p) /* empty */
#define NGH(e,g,h) /* empty */
#define NGP(e,g,p) /* empty */
#define NHP(e,h,p) /* empty */
#define NGHP(e,g,h,p) /* empty */

#define ID(e) /* empty */
#define IS(e) /* empty */
#define ISG(e,g) /* empty */


#endif

