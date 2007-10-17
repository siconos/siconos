/* Siconos-Kernel version 2.1.1, Copyright INRIA 2005-2007.
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

/*! \file

Typedef for simulation-related objects
*/

#ifndef SimulationTypedef_H
#define SimulationTypedef_H

#include <vector>
#include <map>
#include <set>

// ================== Objects to handle UnitaryRelations ==================

#include "UnitaryRelationsSet.h"

/** Map of SiconosMatrices with a UnitaryRelations as a key - Used for diagonal block-terms in assembled matrices of LCP etc ...*/
typedef std::map< UnitaryRelation* , SiconosMatrix*>  MapOfUnitaryMatrices;

/** Iterator through a MapOfUnitaryMatrices */
typedef MapOfUnitaryMatrices::iterator UnitaryMatrixColumnIterator ;

/** Const iterator through a MapOfUnitaryMatrices */
typedef MapOfUnitaryMatrices::const_iterator ConstUnitaryMatrixColumnIterator;

/** Map of MapOfUnitaryMatrices with a UnitaryRelation as a key - Used for extra-diagonal block-terms in assembled matrices of LCP etc ..*/
typedef std::map< UnitaryRelation* , MapOfUnitaryMatrices >  MapOfMapOfUnitaryMatrices;

/** Iterator through a MapOfMapOfUnitaryMatrices */
typedef MapOfMapOfUnitaryMatrices::iterator UnitaryMatrixRowIterator ;

/** Const iterator through a MapOfMapOfUnitaryMatrices */
typedef MapOfMapOfUnitaryMatrices::const_iterator ConstUnitaryMatrixRowIterator ;

/** Map of map of bools, with UnitaryRelations as keys */
typedef std::map< UnitaryRelation* , std::map<UnitaryRelation*, bool> >  MapOfMapOfBool;

/** Vector that contains a sequel of sets of UnitaryRelations*/
typedef std::vector< UnitaryRelationsSet* > VectorOfSetOfUnitaryRelations;

// ================== Objects to handle OSI or DS ==================

#include "OneStepIntegrator.h"

/** Map of SiconosMatrix; key = the related DS*/
typedef std::map<DynamicalSystem*, SiconosMatrix*> MapOfMatrices;

/** Iterator through a map of matrices */
typedef MapOfMatrices::iterator MatIterator;

/** Const iterator through a map of matrices */
typedef MapOfMatrices::const_iterator ConstMatIterator;

/** Map of double; key = the related DS */
typedef std::map<DynamicalSystem*, double> MapOfDouble;

/** Iterator through a map of double */
typedef MapOfDouble::iterator DoubleIterator;

/** Map of bool; key = the related DS */
typedef std::map<DynamicalSystem*, bool> MapOfBool;

/** Vector of OneStepIntegrator */
typedef std::set<OneStepIntegrator*> OSISet;

/** Iterator through vector of OSI*/
typedef OSISet::iterator OSIIterator;

/** Const iterator through vector of OSI*/
typedef OSISet::const_iterator ConstOSIIterator;

/** Return type value for insert function - bool = false if insertion failed. */
typedef std::pair<OSISet::iterator, bool> CheckInsertOSI;

/** A map that links DynamicalSystems and their OneStepIntegrator. */
typedef std::map<DynamicalSystem*, OneStepIntegrator*> DSOSIMap;

/** Iterator through a DSOSIMap. */
typedef DSOSIMap::iterator DSOSIIterator;

/** Const Iterator through a DSOSIMap. */
typedef DSOSIMap::const_iterator DSOSIConstIterator;

// ================== Objects to handle OSNS ==================

//#include "OneStepNSProblem.h"
class OneStepNSProblem;
/** Map of OSNS */
typedef std::map<std::string, OneStepNSProblem* > OneStepNSProblems;

/** Iterator through OneStepNSProblems */
typedef OneStepNSProblems::iterator OSNSIterator;

/** Const iterator through OneStepNSProblems */
typedef OneStepNSProblems::const_iterator ConstOSNSIterator;

// ================== Misc ==================

/** default tolerance value, used to update index sets */
const double DEFAULT_TOLERANCE = 10 * MACHINE_PREC;

#endif
