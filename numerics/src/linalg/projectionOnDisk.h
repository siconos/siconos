/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */

#ifndef ProjectionOnDisk_H
#define ProjectionOnDisk_H

/*!\file projectionOnDisk.h
 * \brief functions to project on a cylinder
 */

#if defined(__cplusplus)
extern "C" {
#endif

/** projectionOnDisk Projection onto the positive Disk of radius R  \f$  \{ r,  0
\sqrt(r_1^2+r_2^2) \geq R \} \f$ \param[in,out] r the vector to be projected \param[in] R the
radius of the cone
*/
void projectionOnDisk(double* r, double R);

/** projectionOnGeneralDisk Projection onto the positive Disk of radius R  \f$  \{ r, 0
\sqrt(r_1^2+r_1^2) \geq R \} \f$ \param[in,out] r the vector to be projected \param[in] R the
radius of the cone \param[in] dim dimension of the cylinder
*/
void projectionOnGeneralDisk(double* r, double R, int dim);

#if defined(__cplusplus)
}
#endif

#endif
