#ifndef LINEARTIEC_H
#define LINEARTIEC_H

#include "LinearEC.h"

/** \class LinearTIEC
 *  \brief Linear Time Invariant Equality Constraint
*  \author SICONOS Development Team - copyright INRIA
 *  \version 0.1
 *  \date 17/01/2005
 *
 *
 */

class LinearTIEC : public LinearEC
{
public:

  /** \fn LinearTIEC(void);
   * \brief default constructor
   */
  LinearTIEC();
  virtual ~LinearTIEC();

  LinearTIEC(EqualityConstraintXML*);

  /** \fn void createEqualityConstraint(LagrangianECXML * ecXML)
   *  \brief allows to create the EqualityConstraint with an xml file, or the needed data
   *  \param LagrangianECXML * : the XML object for this EqualityConstraint
   *  \exception RuntimeException
   */
  void createEqualityConstraint(EqualityConstraintXML * ecXML , int number = -1,
                                SiconosMatrix *G = NULL, vector<DSInputOutput*> *dsioVector = NULL);
};

#endif // LINEARTIEC_H

