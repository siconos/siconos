#ifndef EQUALITYCONSTRAINT_H
#define EQUALITYCONSTRAINT_H


#include "DSInputOutput.h"
#include "EqualityConstraintXML.h"

#include "SiconosSharedLibrary.h"
#include "SiconosMatrix.h"
#include "SiconosConst.h"
#include "check.h"
#include <iostream>

/** \class EqualityConstraint
 *  \brief \todo
*  \author SICONOS Development Team - copyright INRIA
 *  \version 0.1
 *  \date 17/01/2005
 *
 *
 */

class EqualityConstraintXML;
class DSInputOutput;
class SiconosMatrix;
class SiconosSharedLibrary;

class EqualityConstraint
{
public:

  EqualityConstraint();
  EqualityConstraint(EqualityConstraintXML*);

  virtual ~EqualityConstraint();

  /** \fn int getNumber(void);
   *  \brief allows to get the number of the EqualityConstraint
   *  \return the value of number
   */
  inline int getNumber(void)
  {
    return this->number;
  }

  /** \fn string getId(void)
   *  \brief allows to get the id of the EqualityConstraint
   *  \return the value of ths id
   */
  inline std::string  getId(void)
  {
    return this->id;
  }

  /** \fn inline string getType()
   *  \brief allows to get the type of a EqualityConstraint
   *  \return string : the type of the EqualityConstraint
   */
  inline std::string  getType()
  {
    return this->type;
  }

  /** \fn void setNumber(int)
   *  \brief allows to set the value of number
   *  \param int number : an integer to set the value of number
   */
  inline void setNumber(int number)
  {
    this->number = number;
  }

  /** \fn void setId(string)
   *  \brief allows to set the value of id
   *  \param string id : a string to set the value of id
   */
  inline void setId(std::string  id)
  {
    this->id = id;
  }


  /** \fn vector<DSInputOutput*> getDSInputOutputs(void)
   *  \brief allows to get all the DSInputOutput of the EqualityConstraint
   *  \return the vector of DSInputOutput
   */
  std::vector<DSInputOutput*> getDSInputOutputs(void);

  /** \fn DSInputOutput* getDSInputOutput(const int&)
   *  \brief allows to get one specific DSInputOutput, with its place in the vector of DSInputOutput
   *  \param int : the place of the DSInputOutput in the vector of DSInputOutput of the EqualityConstraint
   *  \return DSInputOutput* : dsioVector[ i ] DSInputOutput
   */
  DSInputOutput* getDSInputOutput(const unsigned int&);

  /** \fn void setDSInputOutputs(vector<DSInputOutput*>)
   *  \brief allows to set all the DSInputOutputs of the EqualityConstraint
   *  \param vector<DSInputOutput*> : the vector to set
   */
  void setDSInputOutputs(std::vector<DSInputOutput*>);

  /** \fn void addDSInputOutput(DSInputOutput*)
   *  \brief allows to add the DSInputOutput to the EqualityConstraint
   *  \param DSInputOutput* : the DSInputOutput to add
   */
  void addDSInputOutput(DSInputOutput*);


  /** \fn inline SiconosMatrix* getGPtr()
   *  \brief allows to get G matrix of the EqualityConstraint
   *  \return SiconosMatrix* : the matrix G of the EqualityConstraint
   */
  inline SiconosMatrix* getGPtr()
  {
    return &(this->G);
  }

  /** \fn void setG(SiconosMatrix&)
   *  \brief allows to set the value of G
   *  \param SiconosMatrix& : the matrix to set for G
   */
  inline void setG(SiconosMatrix& G)
  {
    this->G = G;
  }


  /** \fn inline EqualityConstraintXML* getEqualityConstraintXML()
   *  \brief allows to get the EqualityConstraintXML* of the EqualityConstraint
   *  \return a pointer on the EqualityConstraintXML of the EqualityConstraintXML
   */
  inline EqualityConstraintXML* getEqualityConstraintXML()
  {
    return this->ecXML;
  }

  /** \fn inline void setEqualityConstraintXML(EqualityConstraintXML *rxml)
   *  \brief allows to set the EqualityConstraintXML* of the EqualityConstraint
   *  \param EqualityConstraintXML* : the pointer to set
   */
  inline void setEqualityConstraintXML(EqualityConstraintXML *ecxml)
  {
    this->ecXML = ecxml;
  }

  /** \fn void createEqualityConstraint(LagrangianECXML * ecXML)
   *  \brief allows to create the EqualityConstraint with an xml file, or the needed data
   *  \param LagrangianECXML * : the XML object for this EqualityConstraint
   *  \exception RuntimeException
   */
  void createEqualityConstraint(EqualityConstraintXML * ecXML , int number = -1,
                                SiconosMatrix *G = NULL, std::vector<DSInputOutput*> *dsioVector = NULL);

  /** \fn void saveEqualityConstraintToXML()
   *  \brief copy the data of the EqualityConstraint to the XML tree
   */
  void saveEqualityConstraintToXML();

  /** \fn void display()
   *  \brief print the data to the screen
   */
  void display() const;

  //////////////////////////

  /** \fn void computeOutput(double time);
   *  \brief default function to compute y
   *  \param double : current time
   *  \exception RuntimeException
   */
  virtual void computeOutput(double time);

  /** \fn void computeInput(double time);
   *  \brief default function to compute r
   *  \param double : current time
   *  \exception RuntimeException
   */
  virtual void computeInput(double time);

  /** \fn void setComputeOutputFunction(string pluginPath, string functionName)
   *  \brief allow to set a specified function to compute output
   *  \param string : the complete path to the plugin
   *  \param string : the function name to use in this plugin
   *  \exception SiconosSharedLibraryException
   */
  virtual void setComputeOutputFunction(std::string  pluginPath, std::string functionName);

  /** \fn void setComputeInputFunction(string pluginPath, string functionName)
   *  \brief allow to set a specified function to compute output
   *  \param string : the complete path to the plugin
   *  \param string : the function name to use in this plugin
   *  \exception SiconosSharedLibraryException
   */
  virtual void setComputeInputFunction(std::string pluginPath, std::string functionName);

  ///////////////////////


protected :
  /** \fn void fillEqualityConstraintWithEqualityConstraintXML()
   *  \brief uses the EqualityConstraintXML of the EqualityConstraint to fill the fields of this EqualityConstraint
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
  SiconosMatrix G;


  //////////////////////////////
  /** class for manage plugin (open, close librairy...) */
  SiconosSharedLibrary cShared;


  /** \fn void init()
   *  \brief initialise value of a Dynamical System
   */
  virtual void init();


  ////////////////////////
  /* contains the name of the plugin used for computeInput */
  std::string  computeInputName;
  /* contains the name of the plugin used for computeOutput */
  std::string  computeOutputName;

  /** \fn void (*computeOutputPtr)(double* xPtr, double* time, double* lambdaPtr, double* yPtr)
   *  \brief computes y
   *  \param double* xPtr : the pointer to the first element of the vector x
   *  \param double* time : the current time
   *  \param double* lambdaPtr : the pointer to the first element of the vector lambda
   *  \param double* yPtr : the pointer to the first element of the vector y (in-out parameter)
   */
  void (*computeOutputPtr)(double* xPtr, double* time, double* lambdaPtr, double* yPtr);

  /** \fn void (*computeInputPtr)(double* xPtr, double* time, double* lambdaPtr, double* rPtr)
   *  \brief computes r
   *  \param double* xPtr : the pointer to the first element of the vector x
   *  \param double* time : the current time
   *  \param double* lambdaPtr : the pointer to the first element of the vector lambda
   *  \param double* rPtr : the pointer to the first element of the vector r (in-out parameter)
   */
  void (*computeInputPtr)(double* xPtr, double* time, double* lambdaPtr, double* rPtr);

private :

};

#endif // EQUALITYCONSTRAINT_H

