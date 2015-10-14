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
// @file kernel/SomeDefinitions.cpp
// @author Pierre-Brice Wieber
//
// Affiliation(s): INRIA, team BIPOP
//
// Email(s): Pierre-Brice.Wieber@inria.fr
//
// @brief Definitions of the global variables : ActiveDOF and BilateralContact
//

#include <iostream>
#include <vector>
#include <string>

#include "KernelSomeDefinitions.hpp"

/**
 * Global variable storing the enable degrees of freedom, if NULL every joints enabled
 *
 * \var ACTIVEDOF (int) <em> dim depends on the user </em>
 */
std::vector<int> ACTIVEDOF(0);

/**
 * Global variable storing if some contacs are bilateral, if NULL all contacts are unilateral
 *
 * \var BILATERALCONTACT (int) <em> dim NCONT</em>
 */
std::vector<int> BILATERALCONTACT(0);


/***********************************************
 *                 Setters                     *
 ***********************************************/

void setActiveDOF(int * activeDOF, int * size)
{
  if (ACTIVEDOF.size() != 0)
  {
    std::cout << "You are reinitializing ActiveDOF...\n" ;
    ACTIVEDOF.resize(0);
  }

  for (int i = 0; i < *size; i++)
    ACTIVEDOF.push_back(activeDOF[i]);
}

void setBilateralContact(int * bilateralContact, int * size)
{
  if (BILATERALCONTACT.size() != 0)
  {
    std::cout << "You are reinitializing BilateralContact...\n" ;
    BILATERALCONTACT.resize(0);
  }

  for (int i = 0; i < *size; i++)
    BILATERALCONTACT.push_back(bilateralContact[i]);
}

/***********************************************
 *                 Getters                     *
 ***********************************************/

void getActiveDOFsize(int * size)
{
  *size = ACTIVEDOF.size();
}

void getActiveDOF(int * activeDOF)
{
  for (unsigned int i = 0; i < ACTIVEDOF.size(); i++)
    activeDOF[i] = ACTIVEDOF[i];
}

void getBilateralContactsize(int * size)
{
  *size = BILATERALCONTACT.size();
}

void getBilateralContact(int * bilateralContact)
{
  for (unsigned int i = 0; i < BILATERALCONTACT.size(); i++)
    bilateralContact[i] = BILATERALCONTACT[i];
}

