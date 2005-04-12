#ifndef LCP_H
#define LCP_H

#include "OneStepNSProblem.h"
//#include "LCPStructure.h"

#include "SiconosMatrix.h"
//#include "SiconosVector.h"
#include "NewSiconosVector.h"
#include "SimpleVector.h"
//#include "LCPXML.h"


//using namespace std;


class OneStepNSProbem;

/** \class LCP
 *  \brief This class is devoted to the formalization and the resolution of the
 * Linear Complementarity Problem (LCP)
*  \author SICONOS Development Team - copyright INRIA
 *  \version 1.0
 *  \date (Creation) Apr 26, 2004
 *
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
