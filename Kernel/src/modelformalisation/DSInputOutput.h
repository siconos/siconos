//$Id: DSInputOutput.h,v 1.8 2005/03/14 16:05:26 jbarbier Exp $
#ifndef DSIO_H
#define DSIO_H

#include "DynamicalSystem.h"
#include "DSInputOutputXML.h"

#include "SiconosConst.h"

using namespace std;

class Interaction;
class DynamicalSystem;
class DSInputOutputXML;

/** \class DSInputOutput
 *  \brief this class contains data for a specific DynamicalSystem about
 *         { y = H(x, t)
 *         { R = G(r)
 *  \author JEAN MICHEL BARBIER
 *  \version 1.0
 *  \date (Creation) Jan 14, 2005
 *
 * $Date: 2005/03/14 16:05:26 $
 * $Revision: 1.8 $
 * $Author: jbarbier $
 * $Source: /CVS/Siconos/SICONOS/src/modelformalisation/DSInputOutput.h,v $
 *
 *
 *   \warning
 */
class DSInputOutput
{
public:

  /** \fn DSInputOutput()
   *  \brief default constructor
   */
  DSInputOutput();

  /** \fn DSInputOutput(DSInputOutputXML*)
   *  \brief constructor with XML object of the DSInputOutput
   *  \param DSInputOutputXML* : the XML object corresponding
   */
  DSInputOutput(DSInputOutputXML*);

  virtual ~DSInputOutput();


  /** \fn int getNumber(void);
   *  \brief allows to get the number of the EqualityConstraint
   *  \return the value of number
   */
  inline int getNumber(void) const
  {
    return this->number;
  }

  /** \fn string getId(void)
   *  \brief allows to get the id of the EqualityConstraint
   *  \return the value of ths id
   */
  inline string getId(void) const
  {
    return this->id;
  }

  /** \fn void setNumber(int)
   *  \brief allows to set the value of number
   *  \param int number : an integer to set the value of number
   */
  inline void setNumber(const int number)
  {
    this->number = number;
  }

  /** \fn void setId(string)
   *  \brief allows to set the value of id
   *  \param string id : a string to set the value of id
   */
  inline void setId(const string id)
  {
    this->id = id;
  }

  /** \fn inline RelationXML* getDSInputOutputXML()
   *  \brief allows to get the DSInputOutputXML* of the DSInputOutput
   *  \return a pointer on the DSInputOutputXML of the DSInputOutput
   */
  inline DSInputOutputXML* getDSInputOutputXML() const
  {
    return this->dsioxml;
  }

  /** \fn inline void setDSInputOutputXML(DSInputOutputXML *rxml)
   *  \brief allows to set the DSInputOutputXML* of the DSInputOutput
   *  \param DSInputOutputXML* : the pointer to set
   */
  inline void setDSInputOutputXML(DSInputOutputXML *rxml)
  {
    this->dsioxml = rxml;
  }

  /** \fn inline string getType()
   *  \brief allows to get the type of the DSInputOutput
   *  \return string : the type of the DSInputOutput
   */
  inline string getType() const
  {
    return this->dsioType;
  }


  /** \fn inline SiconosMatrix* getHPtr()
   *  \brief allows to get H matrix of the DSInputOutput
   *  \return SiconosMatrix* : the matrix H of the DSInputOutput
   */
  inline SiconosMatrix* getHPtr()
  {
    return &(this->H);
  }

  /** \fn void setH(SiconosMatrix&)
   *  \brief allows to set the value of H
   *  \param SiconosMatrix& : the matrix to set for H
   */
  inline void setH(SiconosMatrix& H)
  {
    this->H = H;
  }

  /** \fn void display()
   *  \brief print the data to the screen
   */
  void display() const;

  /** \fn vector<DynamicalSystem*> getDynamicalSystems(void)
   *  \brief allows to get all the DynamicalSystem of the DSInputOutput
   *  \return the vector of DS
   */
  inline vector<DynamicalSystem*> getDynamicalSystems(void) const
  {
    return this->dsVector;
  };

  /** \fn void setDynamicalSystems(vector<DynamicalSystem*>)
   *  \brief allows to set all the DynamicalSystems of the DSInputOutput
   *  \param vector<DynamicalSystem> : the vector to set
   */
  inline void setDynamicalSystems(const vector<DynamicalSystem*> dsVect)
  {
    this->dsVector = dsVect;
  } ;

  ///////////////////////

  /** \fn void saveDSInputOutputToXML()
   *  \brief copy the data of the DSInputOutput to the XML tree
   */
  void saveDSInputOutputToXML();


  /** \fn void createDSInputOutput(DSInputOutputXML * dsioXML)
   *  \brief allows to create the DSInputOutput with an xml file, or the needed data
   *  \param DSInputOutputXML * : the XML object for this DSInputOutput
   *  \exception RuntimeException
   */
  void createDSInputOutput(DSInputOutputXML * dsioXML, int number = -1,
                           SiconosMatrix *H = NULL);


protected:
  /** \fn void fillDSInputOutputWithDSInputOutputXML()
   *  \brief uses the DSInputOutputXML of the DSInputOutput to fill the fields of this DSInputOutput
   *  \exception RuntimeException
   */
  void fillDSInputOutputWithDSInputOutputXML();



  /** the type of the DSInputOutput : LinearDSIO, LagrangianDSIO */
  string dsioType;

  /** this number defines in a single way the DSInputOutput */
  int number;

  /** the name of the DSInputOutput*/
  string id;

  /** the matrix H */
  SiconosMatrix H;

  /** the DSs connected to this DSInputOuput */
  vector<DynamicalSystem*> dsVector;


  /** the object linked this Relation to read XML data */
  DSInputOutputXML *dsioxml;

  //  /** class for manage plugin (open, close librairy...) */
  //  SiconosSharedLibrary cShared;
  //
  //  /* contains the name of the plugin used for computeInput */
  //  string computeInputName;
  //  /* contains the name of the plugin used for computeOutput */
  //  string computeOutputName;
  //
  //  /** \fn void (*computeOutputPtr)(double* xPtr, double* time, double* lambdaPtr, double* yPtr)
  //   *  \brief computes y
  //   *  \param double* xPtr : the pointer to the first element of the vector x
  //   *  \param double* time : the current time
  //   *  \param double* lambdaPtr : the pointer to the first element of the vector lambda
  //   *  \param double* yPtr : the pointer to the first element of the vector y (in-out parameter)
  //   */
  //  void (*computeOutputPtr)(double* xPtr, double* time, double* lambdaPtr, double* yPtr);
  //
  //  /** \fn void (*computeInputPtr)(double* xPtr, double* time, double* lambdaPtr, double* rPtr)
  //   *  \brief computes r
  //   *  \param double* xPtr : the pointer to the first element of the vector x
  //   *  \param double* time : the current time
  //   *  \param double* lambdaPtr : the pointer to the first element of the vector lambda
  //   *  \param double* rPtr : the pointer to the first element of the vector r (in-out parameter)
  //   */
  //  void (*computeInputPtr)(double* xPtr, double* time, double* lambdaPtr, double* rPtr);

  /** \fn void init()
   *  \brief initialise value of a Relation
   */
  virtual void init();

private :


};

#endif // DSIO_H
//$Log: DSInputOutput.h,v $
//Revision 1.8  2005/03/14 16:05:26  jbarbier
//- manual creation of DSInputOutput saving OK
//
//- in progress for EqualityConstraint
//
//Revision 1.7  2005/03/11 15:06:20  jbarbier
//- save to XML methods of EqualityConstraint and DSInputOutput added
//
//- XML loading process modified : Model loads NSDS, then NSDS loads the DynamicalSystems, EqualityConstraints, Interactions; Modle loads Strategy, then Strategy loads TimeDiscretisation, then the Integrators, then the OneStepNSProblem
//
//Revision 1.6  2005/03/10 12:55:19  jbarbier
//- implmentation of the EqualityConstraint and DSInputOutput classes in progress
//    attributes H (DSIO) et G (EC) added in XML and managed in XML objects
//
//Revision 1.5  2005/03/09 15:30:24  jbarbier
//- add of LagrangianEC class
//
//- in progress : implementation of the EqualityConstraint and DSInputOutput - create methods
//
//Revision 1.4  2005/03/08 12:41:36  jbarbier
//- constant variables files modified :
//Some constants added in SiconosConst
//
//all global tag of the modeling tools are in XMLTagsName, other tags are specific to an XML class
//
//Revision 1.3  2005/02/11 17:35:54  charlety
//
//_ little "inspection of code"
//_ basic getters and setters passed inline
//_ getters functions passed const
//
//Revision 1.2  2005/01/25 14:51:46  jbarbier
//- attributes id, type and XML object added to EqualityConstraint
//
//Revision 1.1  2005/01/17 10:56:24  jbarbier
//- classes EqualityConstraint and DSInputOutput added with inherited classes
//
//- classes EqualityConstraintXML and DSInputOutputXML added with inherited classes
//
