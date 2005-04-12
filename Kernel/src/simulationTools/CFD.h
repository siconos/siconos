#ifndef CFD_H
#define CFD_H

#include "OneStepNSProblem.h"
#include "SiconosMatrix.h"
#include "NewSiconosVector.h"
#include "SimpleVector.h"



class OneStepNSProbem;

/** \class CFD
 *  \brief This class is devoted to the formalization and the resolution of the
 * Contact Friction Dual problem (CFD)
*  \author SICONOS Development Team - copyright INRIA
 *  \version 1.0
 *  \date (Creation) Apr 12, 2005
 */
class CFD : public OneStepNSProblem
{
public:
  /** \fn CFD()
   *  \brief default constructor
   */
  CFD();

  /** \fn CFD(OneStepNSProblemXML*)
   *  \brief constructor with XML object of the CFD
   *  \param OneStepNSProblemXML* : the XML object corresponding
   */
  CFD(OneStepNSProblemXML*);

  ~CFD();

  // getter/setter

  /** \fn SimpleVector getW(void)
   *  \brief get vector w of the CFD
   *  \return SimpleVector : value of w
   */
  inline SimpleVector getW(void) const
  {
    return this->w;
  };

  /** \fn SimpleVector getZ(void)
   *  \brief get vector z of the CFD
   *  \return SimpleVector : value of z
   */
  inline SimpleVector getZ(void) const
  {
    return this->z;
  };



  /** \fn int getNLcp(void)
   *  \brief allow to get the size nLcp of the CFD
   *  \return the size nCfd
   */
  inline int getNCfd(void) const
  {
    return this->nCfd;
  };

  /** \fn SiconosMatrix getM(void)
   *  \brief allow to get the SiconosMatrix M of the CFD
   *  \return the SiconosMatrix M
   */
  inline SiconosMatrix getM(void) const
  {
    return this->M;
  };

  /** \fn SimpleVector getQ(void)
   *  \brief get vector q of the CFD
   *  \return SimpleVector : value of q
   */
  inline SimpleVector getQ(void) const
  {
    return this->q;
  };

  /** \fn SiconosMatrix* getMPtr(void)
  *  \brief allow to get the SiconosMatrix* M of the CFD
  *  \return the SiconosMatrix* M
  */
  SiconosMatrix* getMPtr(void);

  /** \fn SimpleVector* getQPtr(void)
   *  \brief get vector q of the CFD
   *  \return SimpleVector* : pointer on q
   */
  SimpleVector* getQPtr(void);


  /** \fn void setNLcp(int)
   *  \brief set the size of the CFD
   *  \param the size
   */
  inline void setNCfd(const int nCfd)
  {
    this->nCfd = nCfd;
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

  /** \fn CFD* convert (OneStepNSProblem* ds)
   *  \brief encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param OneStepNSProblem* : the one step problem which must be converted
   * \return a pointer on the problem if it is of the right type, NULL otherwise
   */
  static CFD* convert(OneStepNSProblem* osnsp);

private:

  /** Size of the CFD */
  int nCfd;

  /** contains the vector w of a CFD system */
  SimpleVector w;

  /** contains the vector z of a CFD system */
  SimpleVector z;

  /** contains the matrix M of a CFD system */
  SiconosMatrix M;

  /** contains the vector q of a CFD system */
  SimpleVector q;
};

#endif // CFD_H
