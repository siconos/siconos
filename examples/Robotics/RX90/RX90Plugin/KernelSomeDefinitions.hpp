//
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
// @file kernel/SomeDefinitions.hpp
// @author Pierre-Brice Wieber
//
// Affiliation(s): INRIA, team BIPOP
//
// Email(s): Pierre-Brice.Wieber@inria.fr
//
// @brief declarations of the setters and getters
//

#ifndef __Kernel_SomeDefinitions_hpp
#define __Kernel_SomeDefinitions_hpp

#ifdef WINDOWS
#define extern __declspec (dllexport)
#endif

/**
 * set the active degrees of freedom
 *
 * @param[in] size (int)
 * @param[in] activeDOF (int) <em> dim size </em>
 */
extern
#ifdef __cplusplus
"C"
#endif
void setActiveDOF(int *activeDOF, int *size);

/**
 * set the contacts that are bilateral
 *
 * @param[in] size (int)
 * @param[in] bilateral(int) <em> dim size </em>
 */
extern
#ifdef __cplusplus
"C"
#endif
void setBilateralContact(int *bilateralContact, int *size);


/**
 * get the size of ACTIVEDOF
 *
 * @param[out] size (int) size of ACTIVEDOF
 */
extern
#ifdef __cplusplus
"C"
#endif
void getActiveDOFsize(int *size);

/**
 * get the active degrees of freedom
 *
 * @param[out] activeDOF (int) <em> dim ActiveDOFsize </em>
 */
extern
#ifdef __cplusplus
"C"
#endif
void getActiveDOF(int *activeDOF);

/**
 * get the size of BILATERALCONTACT
 *
 * @param[out] size (int) size of BILATERALCONTACT
 */
extern
#ifdef __cplusplus
"C"
#endif
void getBilateralContactsize(int *size);

/**
 * get the contacts that are bilateral
 *
 * @param[out] bilateral(int) <em> dim NCONT </em>
 */
extern
#ifdef __cplusplus
"C"
#endif
void getBilateralContact(int *bilateralContact);


#endif /* __Kernel_SomeDefinitions_hpp */
