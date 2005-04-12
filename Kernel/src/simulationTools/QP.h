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
//$Log: QP.h,v $
//Revision 1.26  2005/02/14 09:52:21  charlety
//_ getters / setters put inline
//
//Revision 1.25  2005/02/04 14:52:44  jbarbier
//- Rolling balls in progress (contact is detected)
//
//- time data is given in parameter for computInput and Output in the Relation. Parameter is sent by methods of the OneStepNSProblem
//
//Revision 1.24  2005/01/31 16:26:26  charlety
//
//_ Added a method named "convert" to classes which inherits from another. This is necessary for Python interface, in order to be able to use down-casting mechanism.
//
//Revision 1.23  2005/01/27 13:57:46  jbarbier
//- suppression of old LCP and QP structures
//
//Revision 1.22  2004/09/23 14:09:24  jbarbier
//- modification of the integrators, the attribute r is always optional.
//
//- modification of the LagrangianNonLinearR. computeInput and computeOutput are
//required.
//
//Revision 1.21  2004/09/22 14:11:14  charlety
//
//  _ revision of Doxygen comments in modelstrategy
//
//Revision 1.20  2004/09/10 11:26:17  charlety
//
//_ Integration of the new version of the SiconosVector in the platform. the class simpleVector is used mostly to replace old SiconosVector. When a vector can be composite or simple, like the state of a dynamical system, a pointer on SiconosVector is used, and the vector is initialized simple or composite when the system is initialized.
//
//_ All the tests which worked with the previous version of the vector are OK with the new version.
//
//_ Example SICONOS and bouncingBall are OK
//
//_ some comments have still to be adapted to NewSiconosVector .
//
//_ This version of NewSiconosVector could be called 0.9. some details have to be fixed, it will be done before the end of September.
//
//Revision 1.19  2004/09/09 08:57:44  jbarbier
//- functions createLCP, createQP, createRelay, addMoreau, addAdams, addLsodar,
//createTimeDiscretisation of the Strategy done.
//
//=> all functions to create manually the objects of the platform are done
//
//Revision 1.18  2004/08/12 11:55:18  jbarbier
//- new methods createModel, createNSDS, createStrategy, ...
//they now allow to make the link with upper objects of the platform
//it will be used for the creation of the platform without XML input file
//
//- the createModel method is finished but the attributes of the other objects
//of the platform are missing for the conctruction
//
//Revision 1.17  2004/08/11 14:43:45  jbarbier
//- beginning of the mechanism of creation without XML input file of the objects of the platform with the
//creatObjtect methods
//
//- function saveWToXML for Moreau integrator, and same specific functions to save
//M,q and Q,p for LCP and QP
//
//- function to check coherency of the Model
//
//Revision 1.16  2004/07/29 14:25:40  jbarbier
