//$Id: LCP.h,v 1.35 2005/02/14 09:52:21 charlety Exp $
#ifndef LCP_H
#define LCP_H

#include "OneStepNSProblem.h"
//#include "LCPStructure.h"

#include "SiconosMatrix.h"
//#include "SiconosVector.h"
#include "NewSiconosVector.h"
#include "SimpleVector.h"
//#include "LCPXML.h"


using namespace std;


class OneStepNSProbem;

/** \class LCP
 *  \brief This class is devoted to the formalization and the resolution of the
 * Linear Complementarity Problem (LCP)
 *  \author Jean-Michel Barbier
 *  \version 1.0
 *  \date (Creation) Apr 26, 2004
 *
 * $Date: 2005/02/14 09:52:21 $
 * $Revision: 1.35 $
 * $Author: charlety $
 * $Source: /CVS/Siconos/SICONOS/src/modelstrategy/LCP.h,v $
 *
 * This class is devoted to the formalization and the resolution of the
 * Linear Complementarity Problem (LCP) defined by :
 *  * \f[
 * w =  q + M z
 * \f]
 * \f[
 * w \geq 0, z \geq 0,  z^{T} w =0
 * \f]
 * where
 *    - \f$w \in R^{n} \f$  and \f$z \in R^{n} \f$ are the unknown,
 *    - \f$A \in R^{n \times n } \f$  and \f$q \in R^{n} \f$
 *
 * \todo Correct the computation of A with a correct concatenation process
 * \warning W is directly used in the computation of A instead of W^{-1}
 */
class LCP : public OneStepNSProblem
{
public:
  /** \fn LCP()
   *  \brief default constructor
   */
  LCP();

  /** \fn LCP(OneStepNSProblemXML*)
   *  \brief constructor with XML object of the LCP
   *  \param OneStepNSProblemXML* : the XML object corresponding
   */
  LCP(OneStepNSProblemXML*);

  ~LCP();

  // getter/setter

  /** \fn SimpleVector getW(void)
   *  \brief get vector w of the LCP
   *  \return SimpleVector : value of w
   */
  inline SimpleVector getW(void) const
  {
    return this->w;
  };

  /** \fn SimpleVector getZ(void)
   *  \brief get vector z of the LCP
   *  \return SimpleVector : value of z
   */
  inline SimpleVector getZ(void) const
  {
    return this->z;
  };



  /** \fn int getNLcp(void)
   *  \brief allow to get the size nLcp of the LCP
   *  \return the size nLcp
   */
  inline int getNLcp(void) const
  {
    return this->nLcp;
  };

  /** \fn SiconosMatrix getM(void)
   *  \brief allow to get the SiconosMatrix M of the LCP
   *  \return the SiconosMatrix M
   */
  inline SiconosMatrix getM(void) const
  {
    return this->M;
  };

  /** \fn SimpleVector getQ(void)
   *  \brief get vector q of the LCP
   *  \return SimpleVector : value of q
   */
  inline SimpleVector getQ(void) const
  {
    return this->q;
  };

  /** \fn SiconosMatrix* getMPtr(void)
  *  \brief allow to get the SiconosMatrix* M of the LCP
  *  \return the SiconosMatrix* M
  */
  SiconosMatrix* getMPtr(void);

  /** \fn SimpleVector* getQPtr(void)
   *  \brief get vector q of the LCP
   *  \return SimpleVector* : pointer on q
   */
  SimpleVector* getQPtr(void);


  /** \fn void setNLcp(int)
   *  \brief set the size of the LCP
   *  \param the size
   */
  inline void setNLcp(const int nLcp)
  {
    this->nLcp = nLcp;
  };

  /** \fn void setM(SiconosMatrix*)
   *  \brief allow to set the SiconosMatrix M
   *  \param the SiconosMatrix to set M
   */
  inline void setM(const SiconosMatrix& M)
  {
    this->M = M;
  };

  /** \fn void setq(SimpleVector&)
   *  \brief set vector q
   *  \param SimpleVector& : new value of q
   */
  inline void setQ(const SimpleVector& Q)
  {
    this->q = Q;
  };

  /** \fn void setW(SimpleVector&)
   *  \brief set vector w
   *  \param SimpleVector& : new value of w
   */
  inline void setW(const SimpleVector& W)
  {
    this->w = W;
  };

  /** \fn void setZ(SimpleVector&)
   *  \brief set vector z
   *  \param SimpleVector& : new value of z
   */
  inline void setZ(const SimpleVector& Z)
  {
    this->z = Z;
  };

  /////////////////////////////////

  /** \fn void formalize(void)
   *  \brief Build the matrix M and the vector b from the OneStep integrator (Problem
   * discretized in time) and the set of interactions.
   *  \param double : current time
   *  \return void
   */
  void formalize(double time);


  /** \fn void compute(void)
   *  \brief Compute the unknown z and w and update the Interaction (y and lambda )
   *  \return void
   */
  void compute(void);

  /** \fn void computeM (void)
   *  \brief make the computation of matrix M
   */
  void computeM(void);

  /** \fn void computeQ (void)
   *  \brief make the computation of the SiconosVector q
   *  \param double : current time
   */
  void computeQ(double time);

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

  /** \fn void saveMToXML()
   *  \brief copy the matrix M of the OneStepNSProblem to the XML tree
   *  \exception RuntimeException
   */
  void saveMToXML();

  /** \fn void saveQToXML()
   *  \brief copy the vector q of the OneStepNSProblem to the XML tree
   *  \exception RuntimeException
   */
  void saveQToXML();

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

  /** \fn LCP* convert (OneStepNSProblem* ds)
   *  \brief encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param OneStepNSProblem* : the one step problem which must be converted
   * \return a pointer on the problem if it is of the right type, NULL otherwise
   */
  static LCP* convert(OneStepNSProblem* osnsp);

private:

  /** Size of the LCP */
  int nLcp;

  /** contains the vector w of a LCP system */
  /*SiconosVector*/
  SimpleVector w;

  /** contains the vector z of a LCP system */
  /*SiconosVector*/
  SimpleVector z;

  /** contains the matrix M of a LCP system */
  SiconosMatrix M;

  /** contains the vector q of a LCP system */
  /*SiconosVector*/
  SimpleVector q;

  //  /** contains the data of the LCP, according to siconos/numerics */
  //  LCPStructure LCPMethod;

  //  /** C structure which gives the informations to the solve function */
  //  methode_lcp meth_lcp;

};

#endif // LCP_H
//$Log: LCP.h,v $
//Revision 1.35  2005/02/14 09:52:21  charlety
//_ getters / setters put inline
//
//Revision 1.34  2005/02/04 14:52:44  jbarbier
//- Rolling balls in progress (contact is detected)
//
//- time data is given in parameter for computInput and Output in the Relation. Parameter is sent by methods of the OneStepNSProblem
//
//Revision 1.33  2005/01/31 16:26:25  charlety
//
//_ Added a method named "convert" to classes which inherits from another. This is necessary for Python interface, in order to be able to use down-casting mechanism.
//
//Revision 1.32  2005/01/27 13:57:46  jbarbier
//- suppression of old LCP and QP structures
//
//Revision 1.31  2004/12/08 12:49:38  jbarbier
//- changes in the XML Schema, respect of the recommandations of the W3C
//version 1.1
//
//- changes in all balises DS, Relation, NSLaw, OneStepIntegrator, OneStepNSProblem
//in the XML files into specific names like LagrangianNLDS, LinearSystemDS, ...
//for the DS
//
//Revision 1.30  2004/12/06 10:10:34  jbarbier
//- integration of Numerics and use of Numerics on the bouncing ball sample
//
//- Numerics is now needed to run the bouncing ball sample!
//
//Revision 1.29  2004/09/23 14:09:24  jbarbier
//- modification of the integrators, the attribute r is always optional.
//
//- modification of the LagrangianNonLinearR. computeInput and computeOutput are
//required.
//
//Revision 1.28  2004/09/22 14:11:13  charlety
//
//  _ revision of Doxygen comments in modelstrategy
//
//Revision 1.27  2004/09/10 11:26:16  charlety
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
//Revision 1.26  2004/09/09 08:57:44  jbarbier
//- functions createLCP, createQP, createRelay, addMoreau, addAdams, addLsodar,
//createTimeDiscretisation of the Strategy done.
//
//=> all functions to create manually the objects of the platform are done
//
//Revision 1.25  2004/08/12 11:55:18  jbarbier
//- new methods createModel, createNSDS, createStrategy, ...
//they now allow to make the link with upper objects of the platform
//it will be used for the creation of the platform without XML input file
//
//- the createModel method is finished but the attributes of the other objects
//of the platform are missing for the conctruction
//
//Revision 1.24  2004/08/11 14:43:45  jbarbier
//- beginning of the mechanism of creation without XML input file of the objects of the platform with the
//creatObjtect methods
//
//- function saveWToXML for Moreau integrator, and same specific functions to save
//M,q and Q,p for LCP and QP
//
//- function to check coherency of the Model
//
//Revision 1.23  2004/07/28 07:43:45  jbarbier
//- all methods createObjectOfThePlatform(...) are now existing
//
//Revision 1.22  2004/07/06 14:54:49  acary
//Renaming NSLaw into NonSmoothLaw
//Renaming RelayNSLaw into RelayNSL
//Renaming CCNSLaw into ComplementarityConditionNSL
//Renaming NewtonImpactLaw into NewtonImpactLawNSL
//
//Revision 1.21  2004/06/30 13:35:56  acary
//Formalization and Computation of the LCP with the NewtonImpactLawNSL
//
//Revision 1.20  2004/06/29 15:12:02  acary
//Change in the naming comvention for the LCP
//The LCP Matrix is now denoted by M.
//The LCP Vector is now denoted by q.
//