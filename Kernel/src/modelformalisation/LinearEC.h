#ifndef LINEAREC_H
#define LINEAREC_H

#include "EqualityConstraint.h"

/** \class LinearEC
 *  \brief Linear Equality Constraint
*  \author SICONOS Development Team - copyright INRIA
 *  \version 0.1
 *  \date 17/01/2005
 *
 *
 */

class LinearEC : public EqualityConstraint
{
public:

  /** \fn LinearEC(void);
   * \brief default constructor
   */
  LinearEC();
  virtual ~LinearEC();

  LinearEC(EqualityConstraintXML*);

  /** \fn void createEqualityConstraint(LagrangianECXML * ecXML)
   *  \brief allows to create the EqualityConstraint with an xml file, or the needed data
   *  \param LagrangianECXML * : the XML object for this EqualityConstraint
   *  \exception RuntimeException
   */
  void createEqualityConstraint(EqualityConstraintXML * ecXML , int number = -1,
                                SiconosMatrix *G = NULL, vector<DSInputOutput*> *dsioVector = NULL);
};

#endif // LINEAREC_H

