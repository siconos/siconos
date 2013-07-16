/* Copyright (C) INRIA 1999-2008
**
** This program is free software; you can redistribute it and/or modify it
** under the terms of the GNU General Public License version 2 as published
** by the Free Software Foundation.
**
** This program is distributed in the hope that it will be useful, but
** WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
** Public License for more details.
**
** You should have received a copy of the GNU General Public License along
** with this program; if not, write to the Free Software Foundation, Inc.,
** 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
**
*/
/**
** @file LagrangianModel/SomeDefinitions.h
** @author: Pierre-Brice Wieber
** Affiliation(s): INRIA, team BIPOP
** Email(s): Pierre-Brice.Wieber@inria.fr
**
** @brief Declarations of the setters and getters on the global variables
**
**
*/

#ifndef __LagrangianModel_SomeDefinitions_h
#define __LagrangianModel_SomeDefinitions_h

#ifdef WINDOWS
#define extern __declspec (dllexport)
#endif

/**
 * set the joints max limit
 *
 * @param[in] maxq <em> dim NDOF </em>
 */
SICONOS_EXPORT void setMaxq(double *maxq);

/**
 * get the joints max limit
 *
 * @param[out] maxq (double) <em> dim NDOF </em>
 */
SICONOS_EXPORT void getMaxq(double *maxq);

/**
 * Global variable corresponding to the min joints limit
 *
 * \var MINQ (double) <em> dim NDOF </em>
 */
//extern double * MINQ;

/**
 * set the joints min limit
 *
 * @param[in] minq <em> dim NDOF </em>
 */
SICONOS_EXPORT void setMinq(double *minq);

/**
 * get the joints min limit
 *
 * @param[out] minq (double) <em> dim NDOF </em>
 */
SICONOS_EXPORT void getMinq(double *minq);

/**
 * Global variable storing the name of the lagrangian model used
 *
 * \var LAGRANGIANMODELNAME (string) <em> dim ? </em>
 */
//extern char* LAGRANGIANMODELNAME;

/**
 * set the lagrangian model name
 *
 * @param[in] lagrangianModelName <em> dim </em>
 */
SICONOS_EXPORT void setLagrangianModelName(char *lagrangianModelName);

/**
 * get the name of the lagrangian model
 *
 * @param[out] lagrangianModelName (char) <em> dim </em>
 */
SICONOS_EXPORT void getLagrangianModelName(char *lagrangianModelName);

/**
 * Global variable storing the name of the contact points
 *
 * \var CONTACTNAMES (string) <em> dim NCONT*? </em>
 */
//extern char** CONTACTNAMES;

/**
 * set the name of the contact corresponding to the index given in input
 *
 * @param[in] index (int)
 * @param[out] contactNames
 */
SICONOS_EXPORT void setContactNames(char *contactNames, int *index);

/**
 * get the name of the contact corresponding to the index given in input
 *
 * @param[in] index (int)
 * @param[out] contactNames (string)
 */
SICONOS_EXPORT void getContactNames(char *contactNames, int *index);

/**
 * Global variable storing the solids of contact
 *
 * \var CONTACTNAMES (int) <em> dim NCONT*? </em>
 */
//extern int** CONTACTSOLIDS;

/**
 * set the solids
 *
 * @param[in] contactSolids <em> dim </em>
 */
SICONOS_EXPORT void setContacSolids(int *contactSolids);

/**
 * get the contact solids
 *
 * @param[out] contactSolids (integer) <em> dim NbSolid*2 </em>
 *
 */
SICONOS_EXPORT void getContactSolids(int *contactSolids);

/**
 * get the number of contact solids
 *
 * @param[out] nbContactSolids (integer)
 *
 */
SICONOS_EXPORT void getNbContactSolids(int *nbContactSolids);


//////////////////////////////////////////////////
//            //
// Human 36         //
//            //
//////////////////////////////////////////////////
/**
 * set the model height
 *
 * @param[in] modelSize Model height (in meter) <em> dim 1 </em>
 */
SICONOS_EXPORT void SetModelSize(double *subjectSize);

/**
 * get the model height
 *
 * @param[out] modelSize Model Height (in meter) <em> dim 1 </em>
 */
SICONOS_EXPORT void GetModelSize(double *subjectSize);

/**
 * set the model mass
 *
 * @param[in] modelMass  Model Mass (in kg) <em> dim 1 </em>
 */
SICONOS_EXPORT void SetModelMass(double *subjectMass);

/**
 * get the model mass
 *
 * @param[out] modelMass Model Mass (in kg) <em> dim 1 </em>
 */
SICONOS_EXPORT void GetModelMass(double *subjectMass);

/**
 * set the model anatomical lengths
 *
 * @param[in] anatLengths anatomical lengths (in meter) row vector <em> dim NLANAT </em>
 */
SICONOS_EXPORT void SetAnatomicalLengths(double *anatLengths);

/**
 * get the model anatomical lengths
 *
 * @param[out] anatLengths anatomical lengths (in meter) row vector <em> dim NLANAT </em>
 */
SICONOS_EXPORT void GetAnatomicalLengths(double *anatLengths);

/**
 * set the model tags positions in their attached segment frames
 *
 * @param[in] tag2JointLengths Row vector of the tags positions in their attached segment frames. The positions are given in meter <em> dim (NTAGS - 1)*3 </em>
 */
SICONOS_EXPORT void SetTag2JointLengths(double *tagLengths);

/**
 * get the model tags positions in their attached segment frames
 *
 * @param[out] tag2JointLengths Row vector of the tags positions in their attached segment frames. The positions are given in meter <em> dim (NTAGS - 1)*3 </em>
 */
SICONOS_EXPORT void GetTag2JointLengths(double *tagLengths);


#endif /* __LagrangianModel_SomeDefinitions_h */

