/* Siconos-Kernel, Copyright INRIA 2005-2011.
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
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
*/

#include "TimeDiscretisationXML.hpp"

TimeDiscretisationXML::TimeDiscretisationXML(xmlNodePtr timeDiscretisationNode):
  _rootNode(timeDiscretisationNode), _hNode(NULL), _NNode(NULL), _tkNode(NULL)
{
  xmlNodePtr node;
  if ((node = SiconosDOMTreeTools::findNodeChild(timeDiscretisationNode, TD_H)))
    _hNode = node;

  if ((node = SiconosDOMTreeTools::findNodeChild(timeDiscretisationNode, TD_N)))
    _NNode = node;

  if ((node = SiconosDOMTreeTools::findNodeChild(timeDiscretisationNode, TD_TK)))
    _tkNode = node;

}

TimeDiscretisationXML::TimeDiscretisationXML(const TimeDiscretisationXML & tdXML)
{
  TimeDiscretisationXML(tdXML.getRootNode());
}


void TimeDiscretisationXML::updateTimeDiscretisationXML(xmlNodePtr node, SP::TimeDiscretisation td)
{
  _rootNode = node;
  XMLException::selfThrow("TimeDiscretisationXML::updateTimeDiscretisationXML - \
      The implementation of this method is not finished !!!");
}


