#ifndef RELAY_H
#define RELAY_H

#include "OneStepNSProblem.h"

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

  /** \fn Relay(OneStepNSProblemXML*, Strategy*=NULL)
   *  \brief xml constructor
   *  \param OneStepNSProblemXML* : the XML linked-object
   *  \param Strategy *: the strategy that owns the problem (optional)
   */
  Relay(OneStepNSProblemXML*, Strategy * = NULL);

  ~Relay();

  /** \fn void saveRelationToXML()
   *  \brief copy the data of the OneStepNSProblem to the XML tree
   */
  void saveNSProblemToXML();

  /** \fn Relay* convert (OneStepNSProblem* ds)
   *  \brief encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param OneStepNSProblem* : the one step problem which must be converted
   * \return a pointer on the problem if it is of the right type, NULL otherwise
   */
  static Relay* convert(OneStepNSProblem* osnsp);
};

#endif // RELAY_H
