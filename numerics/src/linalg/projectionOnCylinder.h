/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */

#ifndef ProjectionOnCylinder_H
#define ProjectionOnCylinder_H

/*!\file projectionOnCylinder.h
 * \brief functions to project on a cylinder
 */

#if defined(__cplusplus)
extern "C" {
#endif

/** projectionOnCylinder Projection onto the positive Cylinder of radius R  \f$  \{ r, r_1 \geq
0, 0 \sqrt(r_2^2+r_3^2) \geq R \} \f$ \param[in,out] r the vector to be projected \param[in] R
the radius of the cone
*/
void projectionOnCylinder(double* r, double R);

/** projectionOnGeneralCylinder Projection onto the positive Cylinder of radius R  \f$  \{ r,
r_1 \geq 0, 0 \sqrt(r_2^2+r_3^2) \geq R \} \f$ \param[in,out] r the vector to be projected
\param[in] R the radius of the cone
\param[in] dim dimension of the cylinder
*/
void projectionOnGeneralCylinder(double* r, double R, int dim);

#if defined(__cplusplus)
}
#endif

#endif
