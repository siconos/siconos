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

/*! \file EqualityConstraint.h

*/

#ifndef EQUALITYCONSTRAINT_H
#define EQUALITYCONSTRAINT_H


#include "DSInputOutput.h"
#include "EqualityConstraintXML.h"

#include "SiconosSharedLibrary.h"
#include "SiconosMatrix.h"
#include "SiconosConst.h"
#include "check.h"
#include <iostream>

const std::string LINEAREC = "LinearEC";
const std::string NLINEAREC = "NLinearEC";
const std::string LINEARTIEC = "LinearTIEC";
const std::string LAGRANGIANEC = "LagrangianEC";
const std::string LAGRANGIANLINEAREC = "LagrangianLinearEC";

class EqualityConstraintXML;
class DSInputOutput;
class SiconosMatrix;
class SiconosSharedLibrary;

//! Not fully implemented
/**  \author SICONOS Development Team - copyright INRIA
 *  \version 2.1.1.
 *  \date 17/01/2005
 *
 *
 */
class EqualityConstraint
{
public:

  EqualityConstraint();
  EqualityConstraint(EqualityConstraintXML*);

  virtual ~EqualityConstraint();

  /** allows to get the number of the EqualityConstraint
  *  \return the value of number
  */
  inline int getNumber(void)
  {
    return this->number;
  }

  /** allows to get the id of the EqualityConstraint
  *  \return the value of ths id
  */
  inline std::string  getId(void)
  {
    return this->id;
  }

  /** allows to get the type of a EqualityConstraint
  *  \return string : the type of the EqualityConstraint
  */
  inline std::string  getType()
  {
    return this->type;
  }

  /** allows to set the value of number
  *  \param int number : an integer to set the value of number
  */
  inline void setNumber(int number)
  {
    this->number = number;
  }

  /** allows to set the value of id
  *  \param string id : a string to set the value of id
  */
  inline void setId(std::string  id)
  {
    this->id = id;
  }


  /** allows to get all the DSInputOutput of the EqualityConstraint
  *  \return the vector of DSInputOutput
  */
  std::vector<DSInputOutput*> getDSInputOutputs(void);

  /** allows to get one specific DSInputOutput, with its place in the vector of DSInputOutput
  *  \param int : the place of the DSInputOutput in the vector of DSInputOutput of the EqualityConstraint
  *  \return DSInputOutput* : dsioVector[ i ] DSInputOutput
  */
  DSInputOutput* getDSInputOutput(const unsigned int&);

  /** allows to set all the DSInputOutputs of the EqualityConstraint
  *  \param vector<DSInputOutput*> : the vector to set
  */
  void setDSInputOutputs(std::vector<DSInputOutput*>);

  /** allows to add the DSInputOutput to the EqualityConstraint
  *  \param DSInputOutput* : the DSInputOutput to add
  */
  void addDSInputOutput(DSInputOutput*);


  /** allows to get G matrix of the EqualityConstraint
  *  \return SiconosMatrix* : the matrix G of the EqualityConstraint
  */
  inline SiconosMatrix* getGPtr()
  {
    return G;
  }

  /** allows to set the value of G
  *  \param SiconosMatrix& : the matrix to set for G
  */
  inline void setG(SiconosMatrix& newG)
  {
    *G = newG;
  }


  /** allows to get the EqualityConstraintXML* of the EqualityConstraint
  *  \return a pointer on the EqualityConstraintXML of the EqualityConstraintXML
  */
  inline EqualityConstraintXML* getEqualityConstraintXML()
  {
    return this->ecXML;
  }

  /** allows to set the EqualityConstraintXML* of the EqualityConstraint
  *  \param EqualityConstraintXML* : the pointer to set
  */
  inline void setEqualityConstraintXML(EqualityConstraintXML *ecxml)
  {
    this->ecXML = ecxml;
  }

  /** allows to create the EqualityConstraint with an xml file, or the needed data
  *  \param LagrangianECXML * : the XML object for this EqualityConstraint
  *  \exception RuntimeException
  */
  void createEqualityConstraint(EqualityConstraintXML * ecXML , int number = -1,
                                SiconosMatrix *G = NULL, std::vector<DSInputOutput*> *dsioVector = NULL);

  /** copy the data of the EqualityConstraint to the XML tree
  */
  void saveEqualityConstraintToXML();

  /** print the data to the screen
  */
  void display() const;

  //////////////////////////

  /** default function to compute y
  *  \param double : current time
  *  \exception RuntimeException
  */
  virtual void computeOutput(double time);

  /** default function to compute r
  *  \param double : current time
  *  \exception RuntimeException
  */
  virtual void computeInput(double time);

  /** allow to set a specified function to compute output
  *  \param string : the complete path to the plugin
  *  \param string : the function name to use in this plugin
  *  \exception SiconosSharedLibraryException
  */
  virtual void setComputeOutputFunction(std::string  pluginPath, std::string functionName);

  /** allow to set a specified function to compute output
  *  \param string : the complete path to the plugin
  *  \param string : the function name to use in this plugin
  *  \exception SiconosSharedLibraryException
  */
  virtual void setComputeInputFunction(std::string pluginPath, std::string functionName);

  ///////////////////////


protected :
  /** uses the EqualityConstraintXML of the EqualityConstraint to fill the fields of this EqualityConstraint
  *  \exception RuntimeException
  */
  void fillEqualityConstraintWithEqualityConstraintXML();

  /** the type of the EqualityConstraint : LinearEC, LinerTIEC, LagrangianEC */
  std::string  type;

  /** this number defines in a single way the EqualityConstraint */
  int number;

  /** the name of the EqualityConstraint*/
  std::string  id;

  /** contains a link to the DSInputOutput of the DynamicalSystems */
  std::vector<DSInputOutput*> dsioVector;

  /** the XML object of this EqualityContraint */
  EqualityConstraintXML *ecXML;

  /** the matrix containing the constraint informations */
  SiconosMatrix* G;


  //////////////////////////////
  /** class for manage plugin (open, close librairy...) */
  SiconosSharedLibrary cShared;


  /** Initialises value of a Dynamical System
   */
  virtual void init();


  ////////////////////////
  /* contains the name of the plugin used for computeInput */
  std::string  computeInputName;
  /* contains the name of the plugin used for computeOutput */
  std::string  computeOutputName;

  /** computes y
  *  \param double* xPtr : the pointer to the first element of the vector x
  *  \param double* time : the current time
  *  \param double* lambdaPtr : the pointer to the first element of the vector lambda
  *  \param double* yPtr : the pointer to the first element of the vector y (in-out parameter)
  */
  void (*computeOutputPtr)(double* xPtr, double* time, double* lambdaPtr, double* yPtr);

  /** computes r
  *  \param double* xPtr : the pointer to the first element of the vector x
  *  \param double* time : the current time
  *  \param double* lambdaPtr : the pointer to the first element of the vector lambda
  *  \param double* rPtr : the pointer to the first element of the vector r (in-out parameter)
  */
  void (*computeInputPtr)(double* xPtr, double* time, double* lambdaPtr, double* rPtr);

private :

};

#endif // EQUALITYCONSTRAINT_H

