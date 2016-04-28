// Copyright (C) INRIA 1999-2008
//
// This program is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License version 2 as published
// by the Free Software Foundation.
//
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
// Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with this program; if not, write to the Free Software Foundation, Inc.,
// 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
//%
// @file ActuationModel/NoDynamics/SomeDefinitions.hpp
// @author Pierre-Brice Wieber
//
// Affiliation(s): INRIA, team BIPOP
//
// Email(s): Pierre-Brice.Wieber@inria.fr
//
// @brief Declarations of the setters and getters on the global variables
//

#ifndef __ActuationModel_SomeDefinitions_hpp
#define __ActuationModel_SomeDefinitions_hpp

#ifdef WINDOWS
#define extern __declspec (dllexport)
#endif

/**
 * set the name of the trajectory
 *
 * @param[in] trajectoryName (string)
 */
extern
#ifdef __cplusplus
"C"
#endif
void setTrajectoryName(char *trajectoryName);

/**
 * set the name defining which contacts are desired during the trajectory
 *
 * @param[in] trajectoryContacts (string)
 */
extern
#ifdef __cplusplus
"C"
#endif
void setTrajectoryContacts(char *trajectoryContacts);

/**
 * set the key positions of the desired trajectory
 *
 * @param[in] size1 (int) number of degrees of freedom
 * @param[in] size2 (int) number of key positions
 * @param[in] trajectoryKeyPositions (double) <em> dim size1*size2 </em>
 */
extern
#ifdef __cplusplus
"C"
#endif
void setTrajectoryKeyPosition(double * trajectoryKeyPositions, int * size1, int * size2);


/**
 * get the name of the trajectory
 *
 * @param[out] trajectoryName (string)
 */
extern
#ifdef __cplusplus
"C"
#endif
void getTrajectoryName(char *trajectoryName);

/**
 * get the name defining which contacts are desired during the trajectory
 *
 * @param[out] trajectoryContacts (string)
 */
extern
#ifdef __cplusplus
"C"
#endif
void getTrajectoryContacts(char *trajectoryContacts);

/**
 * get the number of degrees of freedom
 *
 * @param[in] size (int) number of degrees of freedom
 */
extern
#ifdef __cplusplus
"C"
#endif
void getNDOF(int *size);

/**
 * get the number of key positions
 *
 * @param[in] size (int) number of key positions
 */
extern
#ifdef __cplusplus
"C"
#endif
void getNumberOfKeyPositions(int *size);

/**
 * get the key positions of the desired trajectory
 *
 * @return trajectoryKeyPositions (double) <em> dim NDOF*number of key positions</em>
 */
extern
#ifdef __cplusplus
"C"
#endif
void getTrajectoryKeyPosition(double *trajectoryKeyPositions);

#endif /* __ActuationModel_SomeDefinitions_hpp */
