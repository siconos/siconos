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
// @file ActuationModel/NoDynamics/SomeDefinitions.cpp
// @author Pierre-Brice Wieber
//
// Affiliation(s): INRIA, team BIPOP
//
// Email(s): Pierre-Brice.Wieber@inria.fr
//
// @brief Definitions of the global variables TrajectoryName, TrajectoryContacts,
// TrajectoryKeyPositions,
//

#include <iostream>
#include <vector>
#include <string>

#include "ActuationModelSomeDefinitions.hpp"

#ifdef boost_compil
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#endif


/**
 * Global variable storing the number of actuators
 *
 * \var NB_ACTUATORS (int)
 */
int NB_ACTUATORS = 0;

/**
 * Global variable storing the name of the desired trajectory
 *
 * \var TRAJECTORYNAME (string)
 */
std::string TRAJECTORYNAME;

/**
 * Global variable storing a name allowing to chose the desired contacts in the trajectory
 *
 * \var TRAJECTORYCONTACTS (string)
 */
std::string TRAJECTORYCONTACTS;


/**
 * Global variable storing the desired key positions defining the trajectory
 *
 * \var TRAJECTORYKEYPOSITIONS (double) <em> dim NDOF*nbPos </em>
 */
#ifdef boost_compil
boost::numeric::ublas::matrix<double> TRAJECTORYKEYPOSITIONS(0, 0);
#else
std::vector< std::vector<double> > TRAJECTORYKEYPOSITIONS(0);
#endif


/***********************************************
 *                 Setters                     *
 ***********************************************/

void setNbActuators(int *nbActuators)
{
  NB_ACTUATORS = *nbActuators;
}

void setTrajectoryName(char * trajectoryName)
{
  if (!TRAJECTORYNAME.empty())
  {
    std::cout << "You are reinitializing TrajectoryName...\n" ;
    TRAJECTORYNAME.resize(strlen(trajectoryName));
  }

  TRAJECTORYNAME.replace(0, TRAJECTORYNAME.length(), trajectoryName);
}

void setTrajectoryContacts(char * trajectoryContacts)
{
  if (!TRAJECTORYCONTACTS.empty())
  {
    std::cout << "You are reinitializing TrajectoryContacts...\n" ;
    TRAJECTORYCONTACTS.resize(strlen(trajectoryContacts));
  }

  TRAJECTORYCONTACTS.replace(0, TRAJECTORYCONTACTS.length(), trajectoryContacts);
}

#ifdef boost_compil
void setTrajectoryKeyPosition(double * trajectoryKeyPositions, int * size1, int * size2)
{
  if (TRAJECTORYKEYPOSITIONS.size1()*TRAJECTORYKEYPOSITIONS.size2() != 0)
  {
    std::cout << "You are reinitializing TrajectoryKeyPositions...\n" ;
    TRAJECTORYKEYPOSITIONS.resize(0, 0);
  }

  TRAJECTORYKEYPOSITIONS.resize(*size1, *size2);
  for (int i = 0; i < *size1 * (*size2); i++)
    TRAJECTORYKEYPOSITIONS(i % (*size1), i / (*size1)) = trajectoryKeyPositions[i];
}
#else
void setTrajectoryKeyPosition(double * trajectoryKeyPositions, int * size1, int * size2)
{
  if (TRAJECTORYKEYPOSITIONS.size() != 0)
  {
    std::cout << "You are reinitializing TrajectoryKeyPositions...\n" ;
    for (unsigned int i = 0; i < TRAJECTORYKEYPOSITIONS.size(); i++)
      TRAJECTORYKEYPOSITIONS[i].resize(0);
    TRAJECTORYKEYPOSITIONS.resize(0);
  }

  TRAJECTORYKEYPOSITIONS.resize(*size1);
  for (int i = 0; i < *size1; i++)
  {
    TRAJECTORYKEYPOSITIONS[i].resize(*size2);
    for (int j = 0; j < *size2; j++)
      TRAJECTORYKEYPOSITIONS[i][j] = trajectoryKeyPositions[i * (*size2) + j];
  }
}
#endif


/***********************************************
 *                 Getters                     *
 ***********************************************/

void getNbActuators(int *nbActuators)
{
  *nbActuators = NB_ACTUATORS ;
}

void getTrajectoryName(char * trajectoryName)
{
  TRAJECTORYNAME.copy(trajectoryName, TRAJECTORYNAME.length());
}

void getTrajectoryContacts(char * trajectoryContacts)
{
  TRAJECTORYCONTACTS.copy(trajectoryContacts, TRAJECTORYCONTACTS.length());
}

#ifdef boost_compil
void getNDOF(int *size)
{
  *size = TRAJECTORYKEYPOSITIONS.size1();
}

void getNumberOfKeyPositions(int * size)
{
  *size = TRAJECTORYKEYPOSITIONS.size2();
}

void getTrajectoryKeyPosition(double * trajectoryKeyPositions)
{
  for (unsigned int i = 0; i < TRAJECTORYKEYPOSITIONS.size1()*TRAJECTORYKEYPOSITIONS.size2(); i++)
    trajectoryKeyPositions[i] = TRAJECTORYKEYPOSITIONS(i % (TRAJECTORYKEYPOSITIONS.size1()), i / (TRAJECTORYKEYPOSITIONS.size1()));
}
#else
void getNDOF(int *size)
{
  *size = TRAJECTORYKEYPOSITIONS.size();
}

void getNumberOfKeyPositions(int * size)
{
  *size = TRAJECTORYKEYPOSITIONS[0].size();
}

void getTrajectoryKeyPosition(double * trajectoryKeyPositions)
{
  for (unsigned int i = 0; i < TRAJECTORYKEYPOSITIONS.size(); i++)
    for (unsigned int j = 0; j < TRAJECTORYKEYPOSITIONS[0].size(); j++)
      trajectoryKeyPositions[i * TRAJECTORYKEYPOSITIONS[0].size() + j] = TRAJECTORYKEYPOSITIONS[i][j];
}
#endif

