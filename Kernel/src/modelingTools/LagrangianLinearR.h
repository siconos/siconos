#ifndef LAGRANGIANLINEARRELATION_H
#define LAGRANGIANLINEARRELATION_H

#include "Relation.h"
#include "LagrangianLinearRXML.h"
#include "SimpleVector.h"
#include "CompositeVector.h"

/** \class LagrangianLinearR
 *  \brief Lagrangian Linear Relation
*  \author SICONOS Development Team - copyright INRIA
 *  \version 1.0
 *  \date (Creation) Apr 27, 2004
 *
 *
 */
class LagrangianLinearR : public Relation
{
public:

  /** \fn LagrangianLinearR()
   *  \brief Default constructor
   */
  LagrangianLinearR();

  /** \fn LagrangianLinearR(RelationXML*)
   *  \brief constructor with XML object of the parent class Relation
   *  \param RelationXML* : the XML object corresponding
   */
  LagrangianLinearR(RelationXML*);

  /** \fn LagrangianLinearR(SiconosMatrix, SiconosVector)
   *  \brief constructor with in parameters, the data needed to build this Relation
   *  \param a SiconosMatrix to set h
   *  \param e SiconosVector to set b
   */
  LagrangianLinearR(SiconosMatrix, SimpleVector);
  ~LagrangianLinearR();

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

  /** \fn SiconosMatrix getHPtrRelatingToDS(int position)
   *  \brief getter of the SiconosMatrix H relating to one (the one at the given position in the dsVector of the interaction) of the 2 specific Dynamical System of of the Interaction
   *  \return SiconosMatrix : the H for a DS
   */
  SiconosMatrix getHRelatingToDS(int position);

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
  void saveRelationToXML();

  /** \fn void display()
   *  \brief print the data to the screen
   */
  void display() const;

  /** \fn void createRelation(LagrangianLinearRXML * relationXML,
            SiconosMatrix* H, SiconosVector* b,
            Interaction * interaction)
   *  \brief allows to create the Relation with an xml file, or the needed data
   *  \param LagrangianLinearRXML * : the XML object for this Relation
   *  \param SiconosMatrix* : the matrix H of this Relation
   *  \param SiconosVector* : the vector h of this Relation
   *  \exception RuntimeException
   */
  void createRelation(LagrangianLinearRXML * relationXML,
                      SiconosMatrix* H = NULL, SiconosVector* b = NULL);

  /** \fn LagrangianLinearR* convert (Relation *r)
   *  \brief encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param Relation * : the relation which must be converted
   * \return a pointer on the relation if it is of the right type, NULL otherwise
   */
  static LagrangianLinearR* convert(Relation *r);



protected:
  /** \fn void fillRelationWithRelationXML()
   *  \brief uses the RelationXML of the LagrangianLinearR to fill the fields of this Relation
   *  \exception RuntimeException
   */
  void fillRelationWithRelationXML();


private:
  /** a specific matrix to the LagrangianLinearR */
  SiconosMatrix h;

  /** a specific vector to the LagrangianLinearR */
  /*SiconosVector*/
  SimpleVector b;
};

#endif // LAGRANGIANLINEARRELATION_H
