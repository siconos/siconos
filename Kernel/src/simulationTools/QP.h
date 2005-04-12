#ifndef QP_H
#define QP_H

#include "OneStepNSProblem.h"
//#include "QPStructure.h"
#include "SiconosMatrix.h"
//#include "SiconosVector.h"
#include "NewSiconosVector.h"
//#include "QPXML.h"
#include <iostream>
#include <vector>

class interaction;

//using namespace std;

/** \class QP
 *  \brief It's a way to solve NSDS. It's used in mechanics
*  \author SICONOS Development Team - copyright INRIA
 *  \version 1.0
 *  \date (Creation) Apr 26, 2004
 *
 *
 *
 */
class QP : public OneStepNSProblem
{
public:

  /** \fn QP()
   *  \brief default constructor
   */
  QP();

  /** \fn QP(OneStepNSProblemXML*)
   *  \brief constructor with XML object of the QP
   *  \param OneStepNSProblemXML* : the XML object corresponding
   */
  QP(OneStepNSProblemXML*);

  ~QP();

  // getter/setter
  /** \fn SiconosMatrix getQ(void)
   *  \brief allow to get the SiconosMatrix Q of the QP
   *  \return the SiconosMatrix Q
   */
  inline SiconosMatrix getQ(void) const
  {
    return this->Q;
  };

  /** \fn SimpleVector getP(void)
   *  \brief get vector p of the QP
   *  \return SimpleVector : value of p
   */
  inline SimpleVector getP(void) const
  {
    return this->p;
  };

  /** \fn SiconosMatrix* getQPtr(void)
  *  \brief allow to get the SiconosMatrix* Q of the QP
  *  \return the SiconosMatrix* Q
  */
  SiconosMatrix* getQPtr(void);

  /** \fn SiconosVector* getPPtr(void)
   *  \brief allow to get the SiconosVector* p of the QP
   *  \return the SiconosVector* p
   */
  SimpleVector* getPPtr(void);


  /** \fn void setQ(SiconosMatrix)
   *  \brief allow to set the SiconosMatrix Q
   *  \param the SiconosMatrix to set Q
   */
  inline void setQ(const SiconosMatrix& Q)
  {
    this->Q = Q;
  };

  /** \fn void setP(SimpleVector&)
   *  \brief set vector p
   *  \param SimpleVector& : new value of p
   */
  inline void setP(const SimpleVector& P)
  {
    this->p = p;
  };


  /////////////////////////////

  /** \fn formaliseO(double time)
   *  \brief transform the discretised problem in a problem under numerical form
   *  \param double : current time
   */
  void formalise(double time);

  /** \fn compute()
   *  \brief make the computation so solve the NS problem
   */
  void compute(void);

  /** \fn void fillNSProblemWithNSProblemXML()
   *  \brief uses the OneStepNSProblemXML of the OneStepNSProblem to fill the fields of this OneStepNSProblem
   *  \exception RuntimeException
   */
  void fillNSProblemWithNSProblemXML();

  /** \fn void saveRelationToXML()
   *  \brief copy the data of the OneStepNSProblem to the XML tree
   *  \exception RuntimeException
   */
  void saveNSProblemToXML();

  /** \fn void saveQToXML()
   *  \brief copy the matrix Q of the OneStepNSProblem to the XML tree
   *  \exception RuntimeException
   */
  void saveQToXML();

  /** \fn void savePToXML()
   *  \brief copy the vector p of the OneStepNSProblem to the XML tree
   *  \exception RuntimeException
   */
  void savePToXML();

  /** \fn void display()
   *  \brief print the data to the screen
   */
  void display() const;

  /** \fn void createOneStepNSProblem( OneStepNSProblemXML * osiXML, Strategy * strategy )
   *  \brief allows to create the OneStepNSProblem LCP with an xml file, or the needed data
   *  \param OneStepNSProblemXML * : the XML object for this OneStepNSProblem
   *  \param Strategy * : The NSDS which contains this OneStepNSProblem
   *  \exception RuntimeException
   */
  void createOneStepNSProblem(OneStepNSProblemXML * osiXML, Strategy * strategy = NULL);

  /** \fn QP* convert (OneStepNSProblem* ds)
   *  \brief encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param OneStepNSProblem* : the one step problem which must be converted
   * \return a pointer on the problem if it is of the right type, NULL otherwise
   */
  static QP* convert(OneStepNSProblem* osnsp);

private:
  /** contains the Q matrix of a QP problem */
  SiconosMatrix Q;

  /** contains the p vector of a QP problem */
  SimpleVector p;

  //  /** contains the data of the QP, according to siconos/numerics */
  //  QPStructure QPMethod;

};

#endif // QP_H
