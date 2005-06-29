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
 * y= C x + D \lambda + a
 * \f]
 * and for the output \f$ r\f$ defined by
 * \f[
 * r= B \lambda
 * \f]
 */
class LinearTIR : public Relation
{
public:

  /** \fn LinearTIR();
   *  \brief Basic constructor
   */
  LinearTIR();

  /** \fn LinearTIR(RelationXML*)
   *  \brief constructor with XML object of the LinearTIR
   *  \param RelationXML* : the XML object corresponding
   */
  LinearTIR(RelationXML*);

  /** \fn void LinearTIR(SiconosMatrix* C, SiconosMatrix* D,
                         SiconosMatrix* E, SiconosVector* a,
   *  \brief create the Relation from a set of data
   *  \param SiconosMatrix* : the matrix C
   *  \param SiconosMatrix* : the matrix D (optional)
   *  \param SiconosMatrix* : the matrix E (optional)
   *  \param SiconosVector* : the vector a (optional)
   *  \exception RuntimeException
      */
  LinearTIR(SiconosMatrix* C, SiconosMatrix* D = NULL,
            SiconosMatrix* E = NULL, SiconosVector* a = NULL);

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
  inline void setC(const SiconosMatrix& newValue)
  {
    *C = newValue;
  }

  /** \fn void setCPtr(SiconosMatrix* newPtr)
   *  \brief set C to pointer newPtr
   *  \param SiconosMatrix * newPtr
   */
  inline void setCPtr(SiconosMatrix *newPtr)
  {
    if (isCAllocatedIn) delete C;
    C = newPtr;
    isCAllocatedIn = false;
  }

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
  inline void setD(const SiconosMatrix& newValue)
  {
    *D = newValue;
  }

  /** \fn void setDPtr(SiconosMatrix* newPtr)
   *  \brief set D to pointer newPtr
   *  \param SiconosMatrix * newPtr
   */
  inline void setDPtr(SiconosMatrix *newPtr)
  {
    if (isDAllocatedIn) delete D;
    D = newPtr;
    isDAllocatedIn = false;
  }

  // -- E --

  /** \fn  const SiconosMatrix getE() const
   *  \brief get the value of E
   *  \return SiconosMatrix
   */
  inline const SiconosMatrix getE() const
  {
    return *E;
  }

  /** \fn SiconosMatrix* getEPtr() const
   *  \brief get E
   *  \return pointer on a SiconosMatrix
   */
  inline SiconosMatrix* getEPtr() const
  {
    return E;
  }

  /** \fn void setE (const SiconosMatrix& newValue)
   *  \brief set the value of E to newValue
   *  \param SiconosMatrix newValue
   */
  inline void setE(const SiconosMatrix& newValue)
  {
    *E = newValue;
  }

  /** \fn void setEPtr(SiconosMatrix* newPtr)
   *  \brief set E to pointer newPtr
   *  \param SiconosMatrix * newPtr
   */
  inline void setEPtr(SiconosMatrix *newPtr)
  {
    if (isEAllocatedIn) delete E;
    E = newPtr;
    isEAllocatedIn = false;
  }

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
  inline void setA(const SimpleVector& newValue)
  {
    *a = newValue;
  }

  /** \fn void setAPtr(SimpleVector* newPtr)
   *  \brief set A to pointer newPtr
   *  \param SimpleVector * newPtr
   */
  inline void setAPtr(SimpleVector *newPtr)
  {
    if (isAAllocatedIn) delete a;
    a = newPtr;
    isAAllocatedIn = false;
  }

  // --- OTHER FUNCTIONS ---

  /** \fn void computeOutput()
   *  \brief computes y
   */
  void computeOutput();

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
  /** a matrix specific to the LinearTIR \f$ y= C x + D \lambda + a \f$*/
  SiconosMatrix* C;
  /** a matrix specific to the LinearTIR \f$ y= C x + D \lambda + a \f$*/
  SiconosMatrix* D;
  /** a matrix specific to the LinearTIR \f$ r= E \lambda \f$*/
  SiconosMatrix* E;
  /** a vector specific to the LinearTIR \f$ y= C x + D \lambda + a \f$*/
  SimpleVector* a;

  /** Flags to know if pointers have been allocated in constructors or not*/
  bool isCAllocatedIn;
  bool isDAllocatedIn;
  bool isEAllocatedIn;
  bool isAAllocatedIn;


  //  /** the XML object linked to the LinearTIR to read XML data */
  //  LinearTIRXML *ltirelationxml;
};

#endif // LINEARTIRELATION_H
