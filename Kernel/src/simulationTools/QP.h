#ifndef QP_H
#define QP_H

#include "OneStepNSProblem.h"
#include "QPXML.h"

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

  /** \fn QP(OneStepNSProblemXML*, Strategy*=NULL)
   *  \brief xml constructor
   *  \param OneStepNSProblemXML* : the XML linked-object
   *  \param Strategy *: the strategy that owns the problem (optional)
   */
  QP(OneStepNSProblemXML*, Strategy * = NULL);

  ~QP();

  // --- GETTERS/SETTERS ---

  // --- Q ---
  /** \fn  const SiconosMatrix getQ(void) const
   *  \brief get the value of Q
   *  \return SiconosMatrix
   */
  inline const SiconosMatrix getQ() const
  {
    return *Q;
  }

  /** \fn SiconosMatrix* getQPtr(void) const
   *  \brief get Q
   *  \return pointer on a SiconosMatrix
   */
  inline SiconosMatrix* getQPtr() const
  {
    return Q;
  }

  /** \fn void setQ (const SiconosMatrix& newValue)
   *  \brief set the value of Q to newValue
   *  \param SiconosMatrix newValue
   */
  inline void setQ(const SiconosMatrix& newValue)
  {
    *Q = newValue;
  }

  /** \fn void setQPtr(SiconosMatrix* newPtr)
   *  \brief set Q to pointer newPtr
   *  \param SiconosMatrix * newPtr
   */
  inline void setQPtr(SiconosMatrix *newPtr)
  {
    if (isQAllocatedIn) delete Q;
    Q = newPtr;
    isQAllocatedIn = false;
  }

  // --- P ---
  /** \fn  const SimpleVector getP(void) const
   *  \brief get the value of p, the initial state of the DynamicalSystem
   *  \return SimpleVector
   *  \warning: SiconosVector is an abstract class => can not be an lvalue => return SimpleVector
   */
  inline const SimpleVector getP() const
  {
    return *p;
  }

  /** \fn SimpleVector* getPPtr(void) const
   *  \brief get p, the initial state of the DynamicalSystem
   *  \return pointer on a SimpleVector
   */
  inline SimpleVector* getPPtr() const
  {
    return p;
  }

  /** \fn void setP(const SimpleVector& newValue)
   *  \brief set the value of p to newValue
   *  \param SimpleVector newValue
   */
  inline void setP(const SimpleVector& newValue)
  {
    *p = newValue;
  }

  /** \fn void setPPtr(SimpleVector* newPtr)
   *  \brief set p to pointer newPtr
   *  \param SimpleVector * newPtr
   */
  inline void setPPtr(SimpleVector* newPtr)
  {
    if (isPAllocatedIn) delete p;
    p = newPtr;
    isPAllocatedIn = false;
  }

  // --- OTHER FUNCTIONS ---

  /** \fn formaliseO(const double& time)
   *  \brief transform the discretised problem in a problem under numerical form
   *  \param double : current time
   */
  void formalise(const double& time);

  /** \fn compute()
   *  \brief make the computation so solve the NS problem
   */
  void compute();

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

  /** \fn QP* convert (OneStepNSProblem* ds)
   *  \brief encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param OneStepNSProblem* : the one step problem which must be converted
   * \return a pointer on the problem if it is of the right type, NULL otherwise
   */
  static QP* convert(OneStepNSProblem* osnsp);

private:
  /** \fn QP()
   *  \brief default constructor
   */
  QP();

  /** contains the Q matrix of a QP problem */
  SiconosMatrix* Q;

  /** contains the p vector of a QP problem */
  SimpleVector* p;

  //  /** contains the data of the QP, according to siconos/numerics */
  //  QPStructure QPMethod;

  /** Flags to check wheter pointers were allocated in class constructors or not */
  bool isQAllocatedIn;
  bool isPAllocatedIn;
};

#endif // QP_H
