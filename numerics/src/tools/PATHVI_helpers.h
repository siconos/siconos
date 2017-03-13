/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2017 INRIA.
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

/*! \file PATHVI_helpers.h
 * \brief helpers for using the PATHVI solver
 */

#ifndef PATHVI_HELPERS_H
#define PATHVI_HELPERS_H

#include "SiconosConfig.h"

#include <stddef.h>

/** \struct SN_generic_pathvi_env PATHVI_helpers.h
 * Generic structure for the PATHVI solver*/
typedef struct {
  void* problem;  /**< problem*/
  size_t n;       /**< number of variables*/
  size_t m;       /**< number of polyhedral constraints*/
  double* z;      /**< variable */
  double* F;      /**< function value */
  double* lambda; /**< multipliers for the constraints */
} SN_generic_pathvi_env;

#ifdef HAVE_PATHVI

#include "PATHVI_SDK/include/vi_desc.h"
#define PATHVI_INDX_TYPE int


#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  /** Get the variable value
   * \param desc vi description
   * \param z the variable for PATHVI
   * \return if ok 0, otherwise an error code
   */
  int pathvi_get_z(struct vi_desc *desc, double *z);

  /** Set the variable value
   * \param desc vi description
   * \param z the variable for PATHVI
   * \return if ok 0, otherwise an error code
   */
  int pathvi_set_z(struct vi_desc *desc, double *z);

  /** Get the function value
   * \param desc vi description
   * \param F the function value for PATHVI
   * \return if ok 0, otherwise an error code
   */
  int pathvi_get_F(struct vi_desc *desc, double *F);

  /** Set the function value
   * \param desc vi description
   * \param F the function value for PATHVI
   * \return if ok 0, otherwise an error code
   */
  int pathvi_set_F(struct vi_desc *desc, double *F);

  /** Get the multipliers value
   * \param desc vi description
   * \param lambda the multipliers value for PATHVI
   * \return if ok 0, otherwise an error code
   */
  int pathvi_get_lambda(struct vi_desc *desc, double *lambda);

  /** Set the multipliers value
   * \param desc vi description
   * \param lambda the multipliers value for PATHVI
   * \return if ok 0, otherwise an error code
   */
  int pathvi_set_lambda(struct vi_desc *desc, double *lambda);

  /** Get the name associate with a row (currently "r12" the for 12th row)
   * \param desc vi description
   * \param i the row index
   * \param[out] name the name on output
   * \param len the maximum length to be written in name
   * \return if ok 0, otherwise an error code
   */
  int pathvi_get_row_name(struct vi_desc *desc, int i, char *name, int len);

  /** Get the name associate with a column (currently "c12" the for 12th column)
   * \param desc vi description
   * \param j the column index
   * \param[out] name the name on output
   * \param len the maximum length to be written in name
   * \return if ok 0, otherwise an error code
   */
  int pathvi_get_col_name(struct vi_desc *desc, int j, char *name, int len);

  /** print wrapper
   * \param mode the log mode
   * \param buf the string to print
   * \return if ok 0, otherwise an error code
   */
  void pathvi_print(unsigned mode, const char *buf);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif /* HAVE_PATHVI  */

#endif /* PATHVI_HELPERS_H  */
