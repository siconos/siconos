/* Siconos-Kernel version 2.1.0, Copyright INRIA 2005-2006.
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
/*! \file DSInputOutput.h

*/

#ifndef DSIO_H
#define DSIO_H

#include "SiconosSharedLibrary.h"
#include "SiconosMatrix.h"
#include "SiconosConst.h"
#include "check.h"
#include "SiconosVector.h"
#include "SimpleVector.h"
#include "BlockVector.h"

#include "DynamicalSystem.h"
#include "DSInputOutputXML.h"

#include<string>
#include<iostream>

class DynamicalSystem;
class DSInputOutputXML;
class SiconosMatrix;

const std::string  DefaultComputeInput = "DefaultPlugin:computeInput";
const std::string  DefaultComputeOutput = "DefaultPlugin:computeOutput";
const std::string LINEARDSIO = "LinearDSIO";
const std::string NLINEARDSIO = "NLinearDSIO";
const std::string LAGRANGIANDSIO = "LagrangianDSIO";
const std::string LAGRANGIANLINEARDSIO = "LagrangianLinearDSIO";

//! Not fully implemented
/** this class contains data for a specific DynamicalSystem about
 *         { y = H(x, t)
 *         { R = G(lambda)
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 2.1.0.
 *  \date (Creation) Jan 14, 2005
 *
 *
 *
 *   \warning
 */
class DSInputOutput
{
public:

  /** default constructor
  */
  DSInputOutput();

  /** constructor with XML object of the DSInputOutput
  *  \param DSInputOutputXML* : the XML object corresponding
  */
  DSInputOutput(DSInputOutputXML*);

  virtual ~DSInputOutput();


  /** allows to get the number of the DSIO
  *  \return the value of number
  */
  inline int getNumber(void) const
  {
    return this->number;
  }

  /** allows to get the id of the DSIO
  *  \return the value of ths id
  */
  inline std::string  getId(void) const
  {
    return this->id;
  }

  /** allows to set the value of number
  *  \param int number : an integer to set the value of number
  */
  inline void setNumber(const int number)
  {
    this->number = number;
  }

  /** allows to set the value of id
  *  \param string id : a string to set the value of id
  */
  inline void setId(const std::string  id)
  {
    this->id = id;
  }

  /** allows to get the DSInputOutputXML* of the DSInputOutput
  *  \return a pointer on the DSInputOutputXML of the DSInputOutput
  */
  inline DSInputOutputXML* getDSInputOutputXML() const
  {
    return this->dsioxml;
  }

  /** allows to set the DSInputOutputXML* of the DSInputOutput
  *  \param DSInputOutputXML* : the pointer to set
  */
  inline void setDSInputOutputXML(DSInputOutputXML *rxml)
  {
    this->dsioxml = rxml;
  }

  /** allows to get the type of the DSInputOutput
  *  \return string : the type of the DSInputOutput
  */
  inline std::string  getType() const
  {
    return this->dsioType;
  }


  /** allows to get H matrix of the DSInputOutput
  *  \return SiconosMatrix* : the matrix H of the DSInputOutput
  */
  inline SiconosMatrix* getHPtr()
  {
    return H;
  }

  /** allows to set the value of H
  *  \param SiconosMatrix& : the matrix to set for H
  */
  inline void setH(const SiconosMatrix& newH)
  {
    *H = newH;
  }

  /** print the data to the screen
  */
  void display() const;

  /** allows to get all the DynamicalSystem of the DSInputOutput
  *  \return the vector of DS
  */
  inline std::vector<DynamicalSystem*> getDynamicalSystems(void) const
  {
    return this->dsVector;
  };

  /** allows to set all the DynamicalSystems of the DSInputOutput
  *  \param vector<DynamicalSystem> : the vector to set
  */
  inline void setDynamicalSystems(const std::vector<DynamicalSystem*> dsVect)
  {
    this->dsVector = dsVect;
  } ;

  /** allow to set a specified function to compute output
  *  \param string : the complete path to the plugin
  *  \param string : the function name to use in this plugin
  *  \exception SiconosSharedLibraryException
  */
  virtual void setComputeOutputFunction(std::string  pluginPath, std::string  functionName);

  /** allow to set a specified function to compute output
  *  \param string : the complete path to the plugin
  *  \param string : the function name to use in this plugin
  *  \exception SiconosSharedLibraryException
  */
  virtual void setComputeInputFunction(std::string pluginPath, std::string functionName);

  ///////////////////////

  /** copy the data of the DSInputOutput to the XML tree
  */
  void saveDSInputOutputToXML();


  /** allows to create the DSInputOutput with an xml file, or the needed data
  *  \param DSInputOutputXML * : the XML object for this DSInputOutput
  *  \exception RuntimeException
  */
  void createDSInputOutput(DSInputOutputXML * dsioXML, int number = -1, std::string computeInput = DefaultComputeInput, std::string computeOutput = DefaultComputeOutput);

protected:
  /** uses the DSInputOutputXML of the DSInputOutput to fill the fields of this DSInputOutput
  *  \exception RuntimeException
  */
  void fillDSInputOutputWithDSInputOutputXML();



  /** the type of the DSInputOutput : LinearDSIO, LagrangianDSIO */
  std::string dsioType;

  /** this number defines in a single way the DSInputOutput */
  int number;

  /** the name of the DSInputOutput*/
  std::string id;

  /** the matrix H */
  SiconosMatrix* H;

  /** the DSs connected to this DSInputOuput */
  std::vector<DynamicalSystem*> dsVector;

  /** the object linked this Relation to read XML data */
  DSInputOutputXML *dsioxml;

  /** class for manage plugin (open, close librairy...) */
  SiconosSharedLibrary cShared;

  /* contains the name of the plugin used for computeInput */
  std::string computeInputName;
  /* contains the name of the plugin used for computeOutput */
  std::string computeOutputName;

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

  /** Initialises value of a Relation
   */
  virtual void init();

private :


};

#endif // DSIO_H
