#ifndef LAGRANGIANLINEAREC_H
#define LAGRANGIANLINEAREC_H

#include "EqualityConstraint.h"
#include "LagrangianLinearECXML.h"

/** \class LagrangianLinearEC
 *  \brief Lagrangian EqualityConstraint
*  \author SICONOS Development Team - copyright INRIA
 *  \version 0.1
 *  \date 17/01/2005
 *
 *
 */
class LagrangianLinearEC : public EqualityConstraint
{
public:

  /** \fn LagrangianLinearEC(void);
   * \brief default constructor
   */
  LagrangianLinearEC();
  ~LagrangianLinearEC();

  LagrangianLinearEC(EqualityConstraintXML*);


  /** \fn void createDSInputOutput(EqualityConstraintXML * ecXML)
   *  \brief allows to create the EqualityConstraint with an xml file, or the needed data
   *  \param LagrangianLinearECXML * : the XML object for this EqualityConstraint
   *  \exception RuntimeException
   */
  void createEqualityConstraint(EqualityConstraintXML * ecXML , int number = -1,
                                SiconosMatrix *G = NULL, vector<DSInputOutput*> *dsioVector = NULL);

private:

};

#endif //LAGRANGIANLINEAREC_H 
