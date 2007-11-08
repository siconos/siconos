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
#include "SolverXML.h"
#include "Solver.h" // For DEFAULT_XXX constants

using namespace std;

SolverXML::SolverXML(xmlNodePtr solverNode):
  rootNode(solverNode)
{}

SolverXML::~SolverXML()
{}

string SolverXML::getType() const
{
  return SiconosDOMTreeTools::getStringAttributeValue(rootNode, "type");
}

const bool SolverXML::isBlock() const
{
  if (xmlHasProp(rootNode, (xmlChar *)"isBlock"))
    return SiconosDOMTreeTools::getAttributeValue<bool>(rootNode, "isBlock");
  else return false;
}

unsigned int SolverXML::getMaxIter() const
{
  if (xmlHasProp(rootNode, (xmlChar *)"maxIter"))
    return SiconosDOMTreeTools::getAttributeValue<int>(rootNode, "maxIter");
  else return DEFAULT_ITER;
}

double SolverXML::getTolerance() const
{
  if (xmlHasProp(rootNode, (xmlChar *)"tolerance"))
  {
    string type = SiconosDOMTreeTools::getStringAttributeValue(rootNode, "tolerance");
    return atof(type.c_str());
  }
  else return DEFAULT_TOL;
}

unsigned int SolverXML::getVerbose() const
{
  if (xmlHasProp(rootNode, (xmlChar *)"verbose"))
    return SiconosDOMTreeTools::getAttributeValue<int>(rootNode, "verbose");
  else return DEFAULT_VERBOSE;
}

double SolverXML::getSearchDirection() const
{
  if (xmlHasProp(rootNode, (xmlChar *)"searchDirection"))
  {
    string type = SiconosDOMTreeTools::getStringAttributeValue(rootNode, "searchDirection");
    return atof(type.c_str());
  }
  else return DEFAULT_SEARCHDIR;
}

string SolverXML::getNormType() const
{
  if (xmlHasProp(rootNode, (xmlChar *)"normType"))
    return SiconosDOMTreeTools::getStringAttributeValue(rootNode, "normType");
  else return DEFAULT_NORMTYPE;
}

double SolverXML::getRho() const
{
  if (xmlHasProp(rootNode, (xmlChar *)"Rho"))
  {
    string type = SiconosDOMTreeTools::getStringAttributeValue(rootNode, "Rho");
    return atof(type.c_str());
  }
  else return DEFAULT_RHO;
}

