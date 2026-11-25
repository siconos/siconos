/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */

#ifndef NumericsMatrix_TypesH
#define NumericsMatrix_TypesH
/*!\file NM_types.h
  \brief enum related to NumericsMatrix.
*/

// Declared apart from NumericsMatrix.h to facilitate namespace encapsulation of
// Numerics stuff inside the kernel.

/*! Available types of storage for NumericsMatrix */
typedef enum NumericsMatrix_types {
  NM_DENSE,        /**< dense format */
  NM_SPARSE_BLOCK, /**< sparse block format */
  NM_SPARSE,       /**< compressed column format */
  NM_UNKNOWN,      /**< unset. Used in NM_null */
} NM_types;

#endif
