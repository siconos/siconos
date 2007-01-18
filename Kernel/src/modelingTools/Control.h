/* Siconos-Kernel version 2.0.1, Copyright INRIA 2005-2006.
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

/*! \file Control.h

*/

#ifndef CONTROL
#define CONTROL

#include "NonSmoothDynamicalSystem.h"
#include <string>
#include <map>

class NonSmoothDynamicalSystem;

/**  Base-class to handle control terms in Dynamical Systems.
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 2.0.1.
 *  \date (Creation) January 10, 2007
 *
 */
class Control
{
protected:

  /** Type of Control */
  std::string  type;

  /** size of vector u */
  unsigned int uSize;

  /** "control" term */
  SiconosVector *u;

  /** Matrix coefficient of u */
  SiconosMatrix *T;

  /** Flags to know if pointers have been allocated inside constructors or not */
  AllocationFlagsMap isAllocatedIn;

  /** NonSmoothDynamicalSystem that handles the controlled DynamicalSystems */
  NonSmoothDynamicalSystem* nsds;

  /** Dynamical System concerned by the present control */
  DynamicalSystem * ds;

  /** class for plugin management (open, close librairy...) */
  SiconosSharedLibrary cShared;

  /** Parameters list, last argument of plug-in functions. What are those parameters depends on user's choice.
   *  This a list of pointer to SiconosVector. Each one is identified thanks to a key which is the plug-in name.
   * A flag is also added in the isAllocatedIn map to check inside-class memory allocation for this object.*/
  std::map<std::string, SiconosVector*> parametersList;

  /* the name of the plugin used to compute u and T*/
  std::map<std::string, std::string>  pluginNames;

  /** Flag to check if operators are plugged or not (and thus constant)
   * For example isPlugin["jacobianXF"] = false, means that jacobianXF is constant,
   * then computeJacobianXF does not change its value, and not plugged.*/
  std::map<std::string, bool> isPlugin;

  /** Plug-in to compute u(x,t)
   * @param unsigned int sizeOfU : size of vector u
   * @param unsigned int sizeOfX : size of vector x
   * @param double time : current time
   * @param double* x : pointer to the first element of x
   * @param[in,out] double* u : pointer to the first element of u vector (in-out parameter)
   * @param[in,out] double* param   : a vector of user-defined parameters
   */
  void (*computeUPtr)(unsigned int, unsigned int, double, const double*, double*, double*);

  /** Plug-in to compute T(x)
   * @param unsigned int sizeOfU : size of vector u
   * @param unsigned int sizeOfX : size of vector X
   * @param double* x : pointer to the first element of X
   * @param[in,out] double* T: pointer to the first element of T matrix
   * @param[in,out] double* param   : a vector of user-defined parameters
   */
  void (*computeTPtr)(unsigned int, unsigned int, const double*, double*, double*);

  /** init parameter vector corresponding to id to a SiconosVector* of size 1
   *  \param a string, id of the plug-in
   */
  void initParameter(const std::string);

  /** default constructor
   * \param string: control type, default="default" ...
   */
  Control(const std::string = "default");

public:

  // ===== CONSTRUCTORS =====

  /** copy constructor
   *  \param a Control to be copied
   */
  Control(const Control &);

  // ===== DESTRUCTOR =====

  /** destructor
   */
  virtual ~Control();

  // ===== GETTERS/SETTERS =====

  // --- type of Control ---

  /** get the type of the Control
   *  \return a string
   */
  inline const std::string  getType() const
  {
    return type;
  }

  /** set the type of the Control
   *  \param a string
   */
  inline void setType(const std::string newType)
  {
    type = newType;
  }

  // uSize

  /** to get uSize, size of u
   *  \return the value of uSize
   */
  inline const unsigned int getUSize(void) const
  {
    return uSize;
  }

  /** to set the value of uSize
   *  \param an integer to set the value of uSize
   */
  void setUSize(const unsigned int);

  // ---  U ---

  /** get the value of u, control term
   *  \return SimpleVector
   * \warning: SiconosVector is an abstract class => can not be an lvalue => return SimpleVector
   */
  inline const SimpleVector getU() const
  {
    return *u;
  }

  /** get u, the "control" term
   *  \return pointer on a SiconosVector
   */
  inline SiconosVector* getUPtr() const
  {
    return u;
  }

  /** set the value of u to newValue
   *  \param SiconosVector newValue
   */
  void setU(const SiconosVector&);

  /** set u to pointer newPtr
   *  \param SiconosVector * newPtr
   */
  void setUPtr(SiconosVector *);

  // --- T ---

  /** get the value of T
   *  \return SimpleMatrix
   */
  inline const SimpleMatrix getT() const
  {
    return *T;
  }

  /** get T
   *  \return pointer on a SiconosMatrix
   */
  inline SiconosMatrix* getTPtr() const
  {
    return T;
  }

  /** set the value of T to newValue
   *  \param SiconosMatrix newValue
   */
  void setT(const SiconosMatrix&);

  /** set T to pointer newPtr
   *  \param SiconosMatrix * newPtr
   */
  void setTPtr(SiconosMatrix *newPtr);

  // --- NonSmoothDynamicalSystem ---

  /** get the NonSmoothDynamicalSystem containing this Control.
   *  \return a NonSmoothDynamicalSystem*
   */
  inline NonSmoothDynamicalSystem* getNonSmoothDynamicalSystemPtr() const
  {
    return nsds;
  }

  /** set the NonSmoothDynamicalSystem containing this Control.
   *  \param a NonSmoothDynamicalSystem*
   */
  inline void setNonSmoothDynamicalSystemPtr(NonSmoothDynamicalSystem *newNsds)
  {
    nsds = newNsds;
  }

  // --- DynamicalSystem ---

  /** get the DynamicalSystem concerned by this Control.
   *  \return a DynamicalSystem*
   */
  inline DynamicalSystem* getDynamicalSystemPtr() const
  {
    return ds;
  }

  /** set the DynamicalSystem concerned by this Control.
   *  \param a DynamicalSystem*
   */
  inline void setDynamicalSystemPtr(DynamicalSystem *newDs)
  {
    ds = newDs;
  }

  /** get name of function that computes u (if u from plugin)
   *  \return a string
   */
  inline const std::string getComputeUFunctionName() const
  {
    return pluginNames.find("u")->second;
  }

  /** get name of function that computes T (if T from plugin)
   *  \return a string
   */
  inline const std::string getComputeTFunctionName() const
  {
    return pluginNames.find("T")->second;
  }

  // --- setters for functions to compute plugins ---

  /** to set a specified function to compute u
   *  \param string pluginPath : the complete path to the plugin
   *  \param the string functionName : function name to use in this library
   *  \exception SiconosSharedLibraryException
   */
  void setComputeUFunction(const std::string  pluginPath, const std::string  functionName);

  /** to set a specified function to compute T
   *  \param string pluginPath : the complete path to the plugin
   *  \param the string functionName : function name to use in this library
   *  \exception SiconosSharedLibraryException
   */
  void setComputeTFunction(const std::string  pluginPath, const std::string  functionName);

  // -- parametersList --

  /** get the full map of parameters
   *  \return a map<string,SiconosVector*>
   */
  inline std::map<std::string, SiconosVector*> getParameters() const
  {
    return parametersList;
  };

  /** get the vector of parameters corresponding to plug-in function named id
   *  \return a SimpleVector
   */
  inline const SimpleVector getParameter(const std::string id)
  {
    return *(parametersList[id]);
  };

  /** get the pointer to the vector of parameters corresponding to plug-in function named id
   *  \return a pointer on a SiconosVector
   */
  inline SiconosVector* getParameterPtr(const std::string id)
  {
    return parametersList[id];
  };

  /** set the map for parameters
   *  \param a map<string, SiconosVector*>
   */
  void setParameters(const std::map<std::string, SiconosVector*>&);

  /** set vector corresponding to plug-in function named id to newValue
   *  \param a SiconosVector
   *  \param a string
   */
  void setParameter(const SiconosVector&, const std::string);

  /** set vector corresponding to plug-in function named id to newPtr (!! pointer link !!)
   *  \param a pointer to SiconosVector
   *  \param a string
   */
  void setParameterPtr(SiconosVector *, const std::string);

  // --- compute plugin functions ---

  /** Default function to compute u
   * \param double time : current time
   *  \exception RuntimeException
   */
  virtual void computeU(const double);

  /** function to compute u when x is not those of the current object.
   *  \param double time : current time
   *  \param SiconosVector* : pointer to a x value
   *  \exception RuntimeException
   */
  virtual void computeU(const double,  SiconosVector* xx);

  /** Default function to compute T
   *  \exception RuntimeException
   */
  virtual void computeT();

  // --- isPlugin ---

  /** get isPlugin, map of flags to check if operators are plugged or not
   *  \return a map of bool
   */
  inline const std::map<std::string, bool> getIsPlugin() const
  {
    return isPlugin;
  }

  /** return true if "name" is plugged, else false (ie name is constant)
   *  \return a map of bool
   */
  inline const bool isPlugged(const std::string name)
  {
    return isPlugin[name];
  }

  // ===== MISCELLANEOUS ====

  /** print the data of the dynamical system on the standard output
   */
  virtual void display() const;

};

#endif


