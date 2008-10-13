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

/*! \file RelationTypes.hpp
  \brief enum of the available types and subtypes for relations,
  plugin names ...
*/

#ifndef RELATIONTYPES_HPP
#define RELATIONTYPES_HPP
/** Namespace for user-defined types related to relations */
namespace RELATION
{

/** List of possible DynamicalSystems types*/
enum TYPES
{
  /** First Order */
  FirstOrder,
  /** Lagrangian */
  Lagrangian
};

/** List of possible Relations subtypes*/
enum SUBTYPES
{
  /** non linear */
  NonLinearR,
  /** linear */
  LinearR,
  /** Linear and time invariant */
  LinearTIR,
  /** Scleronomous (lagrangian only) */
  ScleronomousR,
  /** Rheonomous (lagrangian only) */
  RheonomousR,
  /** Compliant (lagrangian only) */
  CompliantR,
  /** */
  Type1R
};

/** The list of all possible plugin names */
enum PluginNames
{
  g, h, G0,
  /** FirstOrder non linear*/
  jacobianH0, jacobianH1, jacobianH2,
  jacobianG0, jacobianG1, jacobianG2,
  /** FirstOrderLinear */
  C, D, F, e, B,
  /** LagrangianScleronomous */
  /** LagrangianRheonomous */
  hDot,
  /** LagrangianCompliant */
  G1, G2
};

/** Container used to check if plugin are plugged or not - PluginNames are used
    to access to PluginBool content. */
typedef std::map<PluginNames, bool> PluginBool;

/** Container used to save the list of plug-in names - PluginNames are used
    to access to PluginBool content. */
typedef std::map<PluginNames, std::string> PluginList;

}
#endif
