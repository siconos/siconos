#ifndef LINEARTIRELATION_H
#define LINEARTIRELATION_H

#include "Relation.h"
#include "LinearTIRXML.h"

/** \class LinearTIR
 *  \brief Linear Time Invariant Relation
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 0.1
 *  \date Apr 27, 2004
 *
 *
 *
 * This class defines and computes the Linear Time Invariant Relation defined by,
 * for the input \f$ y \f$,
 * \f[
 * y= C x + Fu + D \lambda + e
 * \f]
 * and for the output \f$ r\f$ defined by
 * \f[
 * r= B \lambda +a
 * \f]
 */
class LinearTIR : public Relation
{
public:

  /** \fn LinearTIR();
   *  \brief Default constructor
   */
  LinearTIR();

  /** \fn LinearTIR(RelationXML*)
   *  \brief xml constructor
   *  \param LinearTIRXML* : the XML object corresponding
   */
  LinearTIR(RelationXML*);

  /** \fn void LinearTIR(const SiconosMatrix& newC, const SiconosMatrix& newB)
   *  \brief create the Relation from a set of data
   *  \param SiconosMatrix : the matrix C
   *  \param SiconosMatrix : the matrix B
   *  \exception RuntimeException
      */
  LinearTIR(const SiconosMatrix& , const SiconosMatrix&);

  /** \fn void LinearTIR(const SiconosMatrix& newC, const SiconosMatrix& newD,
   *                     const SiconosMatrix& newF, const SimpleVector& newE,
   *                     const SiconosMatrix& newB, const SimpleVector& newA)
   *  \brief create the Relation from a set of data
   *  \param SiconosMatrix : C
   *  \param SiconosMatrix : D
   *  \param SiconosMatrix : F
   *  \param SimpleVectorx : e
   *  \param SiconosMatrix : B
   *  \param SimpleVector : a
   *  \exception RuntimeException
   */
  LinearTIR(const SiconosMatrix& , const SiconosMatrix& ,
            const SiconosMatrix& , const SimpleVector& ,
            const SiconosMatrix& , const SimpleVector&);

  ~LinearTIR();

  // GETTERS/SETTERS

  // -- C --

  /** \fn  const SiconosMatrix getC() const
   *  \brief get the value of C
   *  \return SiconosMatrix
   */
  inline const SiconosMatrix getC() const
  {
    return *C;
  }

  /** \fn SiconosMatrix* getCPtr() const
   *  \brief get C
   *  \return pointer on a SiconosMatrix
   */
  inline SiconosMatrix* getCPtr() const
  {
    return C;
  }

  /** \fn void setC (const SiconosMatrix& newValue)
   *  \brief set the value of C to newValue
   *  \param SiconosMatrix newValue
   */
  void setC(const SiconosMatrix&);

  /** \fn void setCPtr(SiconosMatrix* newPtr)
   *  \brief set C to pointer newPtr
   *  \param SiconosMatrix * newPtr
   */
  void setCPtr(SiconosMatrix *);

  // -- D --

  /** \fn  const SiconosMatrix getD() const
   *  \brief get the value of D
   *  \return SiconosMatrix
   */
  inline const SiconosMatrix getD() const
  {
    return *D;
  }

  /** \fn SiconosMatrix* getDPtr() const
   *  \brief get D
   *  \return pointer on a SiconosMatrix
   */
  inline SiconosMatrix* getDPtr() const
  {
    return D;
  }

  /** \fn void setD (const SiconosMatrix& newValue)
   *  \brief set the value of D to newValue
   *  \param SiconosMatrix newValue
   */
  void setD(const SiconosMatrix&);

  /** \fn void setDPtr(SiconosMatrix* newPtr)
   *  \brief set D to pointer newPtr
   *  \param SiconosMatrix * newPtr
   */
  void setDPtr(SiconosMatrix *);

  // -- F --

  /** \fn  const SiconosMatrix getF() const
   *  \brief get the value of F
   *  \return SiconosMatrix
   */
  inline const SiconosMatrix getF() const
  {
    return *F;
  }

  /** \fn SiconosMatrix* getFPtr() const
   *  \brief get F
   *  \return pointer on a SiconosMatrix
   */
  inline SiconosMatrix* getFPtr() const
  {
    return F;
  }

  /** \fn void setF (const SiconosMatrix& newValue)
   *  \brief set the value of F to newValue
   *  \param SiconosMatrix newValue
   */
  void setF(const SiconosMatrix&);

  /** \fn void setFPtr(SiconosMatrix* newPtr)
   *  \brief set F to pointer newPtr
   *  \param SiconosMatrix * newPtr
   */
  void setFPtr(SiconosMatrix *) ;

  // -- e --

  /** \fn  const SimpleVector getE() const
   *  \brief get the value of e
   *  \return SimpleVector
   */
  inline const SimpleVector getE() const
  {
    return *e;
  }

  /** \fn SimpleVector* getEPtr() const
   *  \brief get e
   *  \return pointer on a SimpleVector
   */
  inline SimpleVector* getEPtr() const
  {
    return e;
  }

  /** \fn void setE (const SimpleVector& newValue)
   *  \brief set the value of e to newValue
   *  \param SimpleVector newValue
   */
  void setE(const SimpleVector&);

  /** \fn void setEPtr(SimpleVector* newPtr)
   *  \brief set E to pointer newPtr
   *  \param SimpleVector * newPtr
   */
  void setEPtr(SimpleVector *);

  // -- B --

  /** \fn  const SiconosMatrix getB() const
   *  \brief get the value of B
   *  \return SiconosMatrix
   */
  inline const SiconosMatrix getB() const
  {
    return *B;
  }

  /** \fn SiconosMatrix* getBPtr() const
   *  \brief get B
   *  \return pointer on a SiconosMatrix
   */
  inline SiconosMatrix* getBPtr() const
  {
    return B;
  }

  /** \fn void setB (const SiconosMatrix& newValue)
   *  \brief set the value of B to newValue
   *  \param SiconosMatrix newValue
   */
  void setB(const SiconosMatrix&);

  /** \fn void setBPtr(SiconosMatrix* newPtr)
   *  \brief set B to pointer newPtr
   *  \param SiconosMatrix * newPtr
   */
  void setBPtr(SiconosMatrix *) ;

  // -- a --

  /** \fn  const SimpleVector getA() const
   *  \brief get the value of a
   *  \return SimpleVector
   */
  inline const SimpleVector getA() const
  {
    return *a;
  }

  /** \fn SimpleVector* getAPtr() const
   *  \brief get a
   *  \return pointer on a SimpleVector
   */
  inline SimpleVector* getAPtr() const
  {
    return a;
  }

  /** \fn void setA (const SimpleVector& newValue)
   *  \brief set the value of a to newValue
   *  \param SimpleVector newValue
   */
  void setA(const SimpleVector&);

  /** \fn void setAPtr(SimpleVector* newPtr)
   *  \brief set a to pointer newPtr
   *  \param SimpleVector * newPtr
   */
  void setAPtr(SimpleVector *);

  // --- OTHER FUNCTIONS ---

  /** \fn void computeOutput()
   *  \brief computes y
   */
  void computeOutput();

  /** \fn void computeInput(double time);
   *  \brief default function to compute lambda
   *  \param double : current time
   *  \exception RuntimeException
   */
  void computeInput();

  /** \fn void saveRelationToXML()
   *  \brief copy the data of the Relation to the XML tree
   */
  void saveRelationToXML();

  /** \fn void display()
   *  \brief print the data to the screen
   */
  void display() const;

  /** \fn LinearTIR* convert (Relation *r)
   *  \brief encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param Relation * : the relation which must be converted
   * \return a pointer on the relation if it is of the right type, NULL otherwise
   */
  static LinearTIR* convert(Relation *r);

private:
  /** Relation is given by: \f$ y= C x + D \lambda + Fu + e \f$*/
  /** and \f$ r = B\lambda\f$ */
  /** */
  SiconosMatrix* C;
  /** */
  SiconosMatrix* D;
  /** */
  SiconosMatrix* F;
  /** */
  SimpleVector* e;
  /** */
  SiconosMatrix* B;
  /** */
  SimpleVector* a;
  /** Flags to know if pointers have been allocated in constructors or not*/
  /* the order is the one of members list above (C,D,F,e,B,a)  */
  std::vector<bool> isAllocatedIn;

  /** the XML object linked to the LinearTIR to read XML data */
  // LinearTIRXML * lTIRxml;
};

#endif // LINEARTIRELATION_H
