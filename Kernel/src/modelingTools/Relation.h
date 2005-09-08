#ifndef RELATION_H
#define RELATION_H

#include "SiconosConst.h"
#include "Interaction.h"
#include "RelationXML.h"
#include "DSInputOutput.h"
#include "check.h"

class Interaction;
class RelationXML;
class DSInputOutput;

/** \class Relation
 *  \brief this class represents relation laws (contact, ...) in an interaction between 2 DS;
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 1.0
 *  \date (Creation) Apr 27, 2004
 *
 *
 *
 *   \warning
 */
class Relation
{
public:

  /** \fn Relation(Interaction* =NULL)
   *  \brief default constructor
   *  \param Interaction*: a pointer to the interaction that owns this relation (optional)
   */
  Relation(Interaction* = NULL);

  /** \fn Relation(RelationXML*)
   *  \brief xml constructor
   *  \param RelationXML* : the XML object corresponding
   *  \param Interaction*: a pointer to the interaction that owns this relation (optional)
   */
  Relation(RelationXML*, Interaction* = NULL);

  /** \fn Relation(const Relation&)
   *  \brief copy constructor
   *  \param a relation to copy
   *  \param Interaction*: a pointer to the interaction that owns this relation (optional)
   *  warning: the interaction link is not copied, set a new one!
   */
  Relation(const Relation&, Interaction* = NULL);

  /** \fn ~Relation()
   *  \brief destructor
   */
  virtual ~Relation();

  /** \fn inline RelationXML* getRelationXML()
   *  \brief allows to get the RelationXML* of the Relation
   *  \return a pointer on the RelationXML of the Relation
   */
  inline RelationXML* getRelationXML()
  {
    return relationxml;
  }

  /** \fn inline void setRelationXML(RelationXML *rxml)
   *  \brief allows to set the RelationXML* of the Relation
   *  \param RelationXML* : the pointer to set
   */
  inline void setRelationXML(RelationXML *rxml)
  {
    relationxml = rxml;
  }

  /** \fn Interaction* getInteractionPtr()
   *  \brief allows to get the Interaction which contains this Relation
   *  \return a pointer on an Interaction
   */
  inline Interaction* getInteractionPtr() const
  {
    return interaction;
  }

  /** \fn void setInteractionPtr(Interaction* i)
   *  \brief set the Interaction which contains this Relation
   */
  inline void setInteractionPtr(Interaction* i)
  {
    interaction = i;
  }

  /** \fn inline string getType()
   *  \brief allows to get the type of the Relation
   *  \return string : the type of the Relation
   */
  inline const std::string  getType() const
  {
    return relationType;
  }

  /** \fn inline string getComputeInputName()
   *  \brief get the name of computeInput function
   *  \return a string
   */
  inline const std::string  getComputeInputName() const
  {
    return computeInputName;
  }

  /** \fn inline string getComputeOutputName()
   *  \brief get the name of computeOutput function
   *  \return a string
   */
  inline const std::string  getComputeOutputName() const
  {
    return computeOutputName;
  }

  /** \fn vector<DSInputOutput*> getDSInputOutputs(void)
   *  \brief allows to get all the DSInputOutput of the Relation
   *  \return the vector of DSInputOutput
   */
  std::vector<DSInputOutput*> getDSInputOutputs(void);

  /** \fn DSInputOutput* getDSInputOutput(const int&)
   *  \brief allows to get one specific DSInputOutput, with its place in the vector of DSInputOutput
   *  \param int : the place of the DSInputOutput in the vector of DSInputOutput of the Relation
   *  \return DSInputOutput* : dsioVector[ i ] DSInputOutput
   */
  DSInputOutput* getDSInputOutput(const unsigned int&);

  /** \fn void setDSInputOutputs(vector<DSInputOutput*>)
   *  \brief allows to set all the DSInputOutputs of the Relation
   *  \param vector<DSInputOutput*> : the vector to set
   */
  void setDSInputOutputs(std::vector<DSInputOutput*>);

  /** \fn void addDSInputOutput(DSInputOutput*)
   *  \brief allows to add the DSInputOutput to the Relation
   *  \param DSInputOutput* : the DSInputOutput to add
   */
  void addDSInputOutput(DSInputOutput*);

  //////////////////////////

  /** \fn void computeOutput(double time);
   *  \brief default function to compute y
   *  \param double : current time
   *  \exception RuntimeException
   */
  virtual void computeOutput(const double&);

  /** \fn void computeFreeOutput(double time);
   *  \brief default function to compute y for the free state
   *  \param double : current time
   *  \exception RuntimeException
   */
  virtual void computeFreeOutput(const double&);

  /** \fn void computeInput(double time);
   *  \brief default function to compute r
   *  \param double : current time
   *  \exception RuntimeException
   */
  virtual void computeInput(const double&);

  /** \fn void setComputeOutputFunction(string pluginPath, string functionName)
   *  \brief allow to set a specified function to compute output
   *  \param string : the complete path to the plugin
   *  \param string : the function name to use in this plugin
   *  \exception SiconosSharedLibraryException
   */
  virtual void setComputeOutputFunction(const std::string&, const std::string&);

  /** \fn void setComputeInputFunction(string pluginPath, string functionName)
   *  \brief allow to set a specified function to compute output
   *  \param string : the complete path to the plugin
   *  \param string : the function name to use in this plugin
   *  \exception SiconosSharedLibraryException
   */
  virtual void setComputeInputFunction(const std::string&, const std::string&);

  ///////////////////////

protected:

  /** type of the Relation */
  std::string  relationType;

  /** the Interaction which contains this Relation */
  Interaction *interaction;

  /** the object linked this Relation to read XML data */
  RelationXML *relationxml;

  /** class for manage plugin (open, close librairy...) */
  SiconosSharedLibrary cShared;

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
  /** contains a link to the DSInputOutput of the DynamicalSystems */
  std::vector<DSInputOutput*> dsioVector;
};

#endif // RELATION_H
