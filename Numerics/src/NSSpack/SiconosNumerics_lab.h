/* Siconos-Numerics version 2.1.0, Copyright INRIA 2005-2006.
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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct
{
  char * nom_method;
  int itermax;
  double tol;
  double mu;
  double k_latin;
} methode_cfp;


typedef struct
{
  char *name;
  int itermax;
  double tol;
  double mu;
  double k_latin;
  double *J1;
  int *ddl_i;
  int *ddl_n;
  int *ddl_tt;
  int *ddl_c;
  int dim_i;
  int dim_c;
  int dim_tt;
  int dim_n;
} method_dfc_2D;


typedef struct
{
  char * name;
  int itermax;
  double tol;
  double k_latin;
} method_lcp;


typedef struct
{
  char * nom_method;
  int itermax;
  double tol;
  double *a;
  double *b;
  double k_latin;
} methode_rp;



typedef struct
{
  char * nom_method;
  int itermax;
  double tol;
  double *a;
  double *b;
  double k_latin;
} methode_rd;




typedef union
{
  methode_cfp cfp;
  methode_cfd cfd;
  methode_lcp lcp;
  methode_rp rp;
  methode_rd rd;

} methode;





