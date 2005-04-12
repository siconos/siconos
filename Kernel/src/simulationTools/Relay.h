#ifndef RELAY_H
#define RELAY_H

#include "OneStepNSProblem.h"
#include <iostream>
#include <vector>
#include "SiconosMatrix.h"
//#include "SiconosVector.h"
#include "NewSiconosVector.h"

//using namespace std;

/** \class Relay
 *  \brief It's a way to solve NSDS. It's used in electronics
*  \author SICONOS Development Team - copyright INRIA
 *  \version 1.0
 *  \date (Creation) Apr 26, 2004
 *
 *
 *
 */
class Relay : public OneStepNSProblem
{
public:
  /** \fn Relay()
   *  \brief default constructor
   */
  Relay();

  /** \fn Relay(OneStepNSProblemXML*)
   *  \brief constructor with XML object of the Relay
   *  \param OneStepNSProblemXML* : the XML object corresponding
   */
  Relay(OneStepNSProblemXML*);

  ~Relay();

  /** \fn void fillNSProblemWithNSProblemXML()
   *  \brief uses the OneStepNSProblemXML of the OneStepNSProblem to fill
   the fields of this OneStepNSProblem
   */
  void fillNSProblemWithNSProblemXML();

  /** \fn void saveRelationToXML()
   *  \brief copy the data of the OneStepNSProblem to the XML tree
   */
  void saveNSProblemToXML();

  /** \fn void createOneStepNSProblem( OneStepNSProblemXML * osiXML, Strategy * strategy )
   *  \brief allows to create the OneStepNSProblem LCP with an xml file, or the needed data
   *  \param OneStepNSProblemXML * : the XML object for this OneStepNSProblem
   *  \param Strategy * : The NSDS which contains this OneStepNSProblem
   *  \exception RuntimeException
   */
  void createOneStepNSProblem(OneStepNSProblemXML * osiXML, Strategy * strategy = NULL);

  /** \fn Relay* convert (OneStepNSProblem* ds)
   *  \brief encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param OneStepNSProblem* : the one step problem which must be converted
   * \return a pointer on the problem if it is of the right type, NULL otherwise
   */
  static Relay* convert(OneStepNSProblem* osnsp);
};

#endif // RELAY_H
