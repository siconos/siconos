#ifndef LAGRANGIANDSIO_H
#define LAGRANGIANDSIO_H

#include "EqualityConstraint.h"
#include "LagrangianECXML.h"

/** \class LagrangianEC
 *  \brief Lagrangian EqualityConstraint
*  \author SICONOS Development Team - copyright INRIA
 *  \version 0.1
 *  \date 17/01/2005
 *
 *
 */
class LagrangianEC : public EqualityConstraint
{
public:

  /** \fn LagrangianEC(void);
   * \brief default constructor
   */
  LagrangianEC();
  ~LagrangianEC();

  LagrangianEC(EqualityConstraintXML*);


  /** \fn void createDSInputOutput(EqualityConstraintXML * ecXML)
   *  \brief allows to create the EqualityConstraint with an xml file, or the needed data
   *  \param LagrangianECXML * : the XML object for this EqualityConstraint
   *  \exception RuntimeException
   */
  void createEqualityConstraint(EqualityConstraintXML * ecXML , int number = -1,
                                SiconosMatrix *G = NULL, vector<DSInputOutput*> *dsioVector = NULL);

private:

};

#endif
