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

  /** \fn LagrangianLinearR(SiconosMatrix*, SiconosVector*)
   *  \brief constructor with in parameters, the data needed to build this Relation
   *  \param a SiconosMatrix to set H
   *  \param a SiconosVector to set b
   */
  LagrangianLinearR(SiconosMatrix*, SimpleVector*);
  ~LagrangianLinearR();

  // --- GETTERS/SETTERS
  // -- H --

  /** \fn  const SiconosMatrix getH() const
   *  \brief get the value of H
   *  \return SiconosMatrix
   */
  inline const SiconosMatrix getH() const
  {
    return *H;
  }

  /** \fn SiconosMatrix* getHPtr() const
   *  \brief get H
   *  \return pointer on a SiconosMatrix
   */
  inline SiconosMatrix* getHPtr() const
  {
    return H;
  }

  /** \fn void setH (const SiconosMatrix& newValue)
   *  \brief set the value of H to newValue
   *  \param SiconosMatrix newValue
   */
  inline void setH(const SiconosMatrix& newValue)
  {
    *H = newValue;
  }

  /** \fn void setHPtr(SiconosMatrix* newPtr)
   *  \brief set H to pointer newPtr
   *  \param SiconosMatrix * newPtr
   */
  inline void setHPtr(SiconosMatrix *newPtr)
  {
    if (isHAllocatedIn) delete H;
    H = newPtr;
    isHAllocatedIn = false;
  }
  // -- b --

  /** \fn  const SimpleVector getB() const
   *  \brief get the value of b
   *  \return SimpleVector
   */
  inline const SimpleVector getB() const
  {
    return *b;
  }

  /** \fn SimpleVector* getBPtr() const
   *  \brief get b
   *  \return pointer on a SimpleVector
   */
  inline SimpleVector* getBPtr() const
  {
    return b;
  }

  /** \fn void setB (const SimpleVector& newValue)
   *  \brief set the value of b to newValue
   *  \param SimpleVector newValue
   */
  inline void setB(const SimpleVector& newValue)
  {
    *b = newValue;
  }

  /** \fn void setBPtr(SimpleVector* newPtr)
   *  \brief set B to pointer newPtr
   *  \param SimpleVector * newPtr
   */
  inline void setBPtr(SimpleVector *newPtr)
  {
    if (isBAllocatedIn) delete b;
    b = newPtr;
    isBAllocatedIn = false;
  }

  /** \fn SiconosMatrix getHPtrRelatingToDS(int position)
   *  \brief getter of the SiconosMatrix H relating to one (the one at the given position in the dsVector of the interaction) of the 2 specific Dynamical System of of the Interaction
   *  \return SiconosMatrix : the H for a DS
   */
  SiconosMatrix getHRelatingToDS(const int&);

  // --- OTHER FUNCTIONS ---

  /** \fn void computeFreeOutput(double time);
   *  \brief default function to compute y for the free state
   *  \param double : current time
   *  \exception RuntimeException
   */
  void computeFreeOutput(const double& time);

  /** \fn void computeOutput(double time);
   *  \brief default function to compute y
   *  \param double : current time
   *  \exception RuntimeException
   */
  void computeOutput(const double& time);

  /** \fn void computeInput(double time);
   *  \brief default function to compute lambda
   *  \param double : current time
   *  \exception RuntimeException
   */
  void computeInput(const double& time);

  /** \fn void saveRelationToXML()
   *  \brief copy the data of the Relation to the XML tree
   *  \exception RuntimeException
   */
  void saveRelationToXML();

  /** \fn void display()
   *  \brief print the data to the screen
   */
  void display() const;

  /** \fn LagrangianLinearR* convert (Relation *r)
   *  \brief encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param Relation * : the relation which must be converted
   * \return a pointer on the relation if it is of the right type, NULL otherwise
   */
  static LagrangianLinearR* convert(Relation *r);

private:
  /** a specific matrix to the LagrangianLinearR */
  SiconosMatrix* H;

  /** a specific vector to the LagrangianLinearR */
  SimpleVector* b;

  /* Flags to control where pointers have been allocated (in or out constructors)*/
  bool isHAllocatedIn;
  bool isBAllocatedIn;


};

#endif // LAGRANGIANLINEARRELATION_H
