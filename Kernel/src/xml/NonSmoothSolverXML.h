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
/*! \file NonSmoothSolverXML
  \brief XML tools for NonSmoothSolver parameters reading/writing
*/

#ifndef __NONSMOOTHSOLVERXML__
#define __NONSMOOTHSOLVERXML__

#include "SiconosDOMTreeTools.h"

/** XML management for NonSmoothSolver
 *
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 3.0.0.
 *   \date 20/12/2005
 *
 * This class is used to read parameters for a NonSmoothSolver in a XML file. \n
 * The parameters are: \n
 *     - the name of the solver/algorithm \n
 *     - a list (vector) of integer parameters \n
 *     - a list (vector) of double parameters \n
 *
 *  In the XML file, NonSmoothSolver must look like:
 *
 \code
   <NonSmoothSolver Name="Lemke">
   <iparam  vectorSize='3'>1 3 4</iparam>
   <dparam  vectorSize='3'>1.3 1e-17 1e-12</dparam>
   </NonSmoothSolver>
 \endcode
 *
 * The tag NonSmoothSolver is optional and, if omitted, solver type and parameters are set during the call to the Numerics driver
 * by reading a default parameters file.
 */
class NonSmoothSolverXML
{
private:

  /** root node named "NonSmoothSolver" - Child of OneStepNSProblem node */
  xmlNodePtr rootNode;

  /** Node for int parameters */
  xmlNodePtr iparamNode;

  /** Node for double parameters */
  xmlNodePtr dparamNode;

  /** Default constructor (private) */
  NonSmoothSolverXML();

  /** Copy constructor (private)*/
  NonSmoothSolverXML(const NonSmoothSolverXML&);

public:

  /** Constructor from xml node
   *  \param : rootNode
   */
  NonSmoothSolverXML(xmlNodePtr);

  /** Destructor
   */
  ~NonSmoothSolverXML();

  /** To get the node corresponding to the NonSmoothSolver
   *  \return pointer to node named "NonSmoothSolver"
   */
  inline xmlNodePtr getRootNode() const
  {
    return rootNode;
  };

  /** To get the name of the solver/algorithm (attribute of root node)
   *  \return a string
   */
  inline std::string getName() const
  {
    return SiconosDOMTreeTools::getStringAttributeValue(rootNode, "Name");
  };

  /** returns true if iparamNode is defined
   *  \return a bool
   */
  inline bool hasIparam() const
  {
    return (iparamNode);
  }

  /** returns true if dparamNode is defined
   *  \return a bool
   */
  inline bool hasDparam() const
  {
    return (dparamNode);
  }

  /** Gets iparam vector values
      \param[in-out] the vector to be filled
   */
  void getIParam(std::vector<int>&);

  /** Gets dparam vector values
      \param[in-out] the vector to be filled
   */
  void getDParam(std::vector<double>&);
};


#endif
