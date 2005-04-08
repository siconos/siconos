#ifndef LAGRANGIANLINEAREC_H
#define LAGRANGIANLINEAREC_H

#include "EqualityConstraint.h"
#include "LagrangianLinearECXML.h"
#include "SimpleVector.h"
#include "CompositeVector.h"

/** \class LagrangianLinearEC
 *  \brief Lagrangian Linear EqualityConstraint
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 1.0
 *  \date (Creation) Apr 27, 2004
 *
 *
 */
class LagrangianLinearEC : public EqualityConstraint
{
public:

  /** \fn LagrangianLinearEC()
   *  \brief Default constructor
   */
  LagrangianLinearEC();

  /** \fn LagrangianLinearEC(EqualityConstraintXML*)
   *  \brief constructor with XML object of the parent class EqualityConstraint
   *  \param EqualityConstraintXML* : the XML object corresponding
   */
  LagrangianLinearEC(EqualityConstraintXML*);

  /** \fn LagrangianLinearEC(SiconosMatrix, SiconosVector)
   *  \brief constructor with in parameters, the data needed to build this Relation
   *  \param a SiconosMatrix to set h
   *  \param e SiconosVector to set b
   */
  LagrangianLinearEC(SiconosMatrix, SimpleVector);
  ~LagrangianLinearEC();

  /** \fn SiconosMatrix getH(void)
   *  \brief getter of the SiconosMatrix h
   *  \return a pointer on the SiconosMatrix h
   */
  inline SiconosMatrix getH(void) const
  {
    return this->h;
  } ;

  /** \fn SimpleVector getB(void)
   *  \brief getter of the SiconosVector b
   *  \return SimpleVector : value of b
   */
  inline /*SiconosVector*/SimpleVector getB(void) const
  {
    return this->b;
  };

  /** \fn SiconosMatrix* getHPtr(void)
   *  \brief getter of the SiconosMatrix* h
   *  \return a pointer on the SiconosMatrix* h
   */
  SiconosMatrix* getHPtr(void);

  /** \fn SiconosVector* getBPtr(void)
   *  \brief getter of the SiconosVector* b
   *  \return a pointer on the SiconosVector b
   */
  SiconosVector* getBPtr(void);

  /** \fn void setH(SiconosMatrix)
   *  \brief setter on the SiconosMatrix h
   *  \param a SiconosMatrix to set h
   */
  inline void setH(const SiconosMatrix &H)
  {
    this->h = H;
  };

  /** \fn void setH(SimpleVector&)
   *  \brief set the vector b
   *  \param SimpleVector& : new value of b
   */
  inline void setB(/*SiconosVector*/SimpleVector& B)
  {
    this->b = B;
  };


  ////////////////////////////

  /** \fn void computeFreeOutput(double time);
   *  \brief default function to compute y for the free state
   *  \param double : current time
   *  \exception RuntimeException
   */
  void computeFreeOutput(double time);

  /** \fn void computeOutput(double time);
   *  \brief default function to compute y
   *  \param double : current time
   *  \exception RuntimeException
   */
  void computeOutput(double time);

  /** \fn void computeInput(double time);
   *  \brief default function to compute lambda
   *  \param double : current time
   *  \exception RuntimeException
   */
  void computeInput(double time);

  /** \fn void saveRelationToXML()
   *  \brief copy the data of the Relation to the XML tree
   *  \exception RuntimeException
   */
  void saveEqualityConstraintToXML();

  /** \fn void display()
   *  \brief print the data to the screen
   */
  void display() const;

  /** \fn void createEqualityConstraint(EqualityConstraintXML * ecXML,
            SiconosMatrix* H, SiconosVector* b)
   *  \brief allows to create the Relation with an xml file, or the needed data
   *  \param LagrangianLinearECXML * : the XML object for this EqualityConstraint
   *  \param SiconosMatrix* : the matrix H of this EqualityConstraint
   *  \param SiconosVector* : the vector h of this EqualityConstraint
   *  \exception RuntimeException
   */
  void createEqualityConstraint(EqualityConstraintXML * ecXML,
                                SiconosMatrix* H = NULL, SiconosVector* b = NULL);

  /** \fn LagrangianLinearEC* convert (EqualityConstraint *ec)
   *  \brief encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param Relation * : the EqualityConstraint which must be converted
   * \return a pointer on the EqualityConstraint if it is of the right type, NULL otherwise
   */
  static LagrangianLinearEC* convert(EqualityConstraint *ec);



protected:
  /** \fn void fillEqualityConstraintWithEqualityConstraintXML()
   *  \brief uses the EqualityConstraintXML of the LagrangianLinearEC to fill the fields of this EqualityConstraint
   *  \exception RuntimeException
   */
  void fillEqualityConstraintWithEqualityConstraintXML();


private:
  /** a specific matrix to the LagrangianLinearEC */
  SiconosMatrix h;

  /** a specific vector to the LagrangianLinearEC */
  /*SiconosVector*/
  SimpleVector b;
};

#endif // LAGRANGIANLINEAREC_H
