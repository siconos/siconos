/* Siconos-Kernel version 3.0.0, Copyright INRIA 2005-2008.
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
/*! \file UnitaryRelationsSet.hpp
 */

#ifndef UnitaryRelationsSET_H
#define UnitaryRelationsSET_H

#include "SiconosSet.hpp"
#include "UnitaryRelation.hpp"

#include "SiconosPointers.hpp"

/** A set of pointers to interactions, sorted in a growing order according to their address */
typedef SiconosSet<UnitaryRelation, double*> UnitaryRelationsSet;
/** Iterator through a set of UnitaryRelations */
typedef std::set<SP::UnitaryRelation, Cmp<UnitaryRelation, double*> >::iterator UnitaryRelationsIterator;
/** const Iterator through a set of UnitaryRelations */
typedef std::set<SP::UnitaryRelation, Cmp<UnitaryRelation, double*> >::const_iterator ConstUnitaryRelationsIterator;
/** return type value for insert function - bool = false if insertion failed. */
typedef std::pair<UnitaryRelationsIterator, bool> CheckInsertUnitaryRelation;


TYPEDEF_SPTR(UnitaryRelationsSet);
#endif
