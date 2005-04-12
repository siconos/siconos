#ifndef LINEARTIRELATION_H
#define LINEARTIRELATION_H

#include "Relation.h"
//#include "RelationXML.h"
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

  ~LinearTIR();

  /** \fn SiconosMatrix getC(void)
   *  \brief allows to get the matrix h of the Relation
   *  \exception to be defined
   *  \return the SiconosMatrix h of the Relation
   */
  inline SiconosMatrix getC(void) const
  {
    return this->C;
  };

  /** \fn SiconosMatrix getD(void)
   *  \brief allows to get the matrix d of the Relation
   *  \exception to be defined
   *  \return the SiconosMatrix d of the Relation
   */
  inline SiconosMatrix getD(void) const
  {
    return this->D;
  };

  /** \fn SiconosMatrix getE(void)
   *  \brief allows to get the matrix E of the Relation
   *  \exception to be defined
   *  \return the SiconosMatrix E of the Relation
   */
  inline SiconosMatrix getE(void) const
  {
    return this->E;
  };

  /** \fn SimpleVector getA(void)
   *  \brief get vector a of the Relation
   *  \exception to be defined
   *  \return SimpleVector : value of a
   */
  inline /*SiconosVector*/SimpleVector getA(void) const
  {
    return this->a;
  };


  /** \fn SiconosMatrix* getCPtr(void)
   *  \brief allows to get the SiconosMatrix* h of the Relation
   *  \exception to be defined
   *  \return the SiconosMatrix* h of the Relation
   */
  SiconosMatrix* getCPtr(void);

  /** \fn SiconosMatrix* getDPtr(void)
   *  \brief allows to get the SiconosMatrix* d of the Relation
   *  \exception to be defined
   *  \return the SiconosMatrix* d of the Relation
   */
  SiconosMatrix* getDPtr(void);

  /** \fn SiconosMatrix* getEPtr(void)
   *  \brief allows to get the SiconosMatrix* E of the Relation
   *  \exception to be defined
   *  \return the SiconosMatrix* E of the Relation
   */
  SiconosMatrix* getEPtr(void);

  /** \fn SiconosVector* getAPtr(void)
   *  \brief allows to get the SiconosVector* a of the Relation
   *  \exception to be defined
   *  \return the SiconosVector* a of the Relation
   */
  SiconosVector* getAPtr(void);

  /** \fn void setC(SiconosMatrix)
   *  \brief allows to set the matrix C of the Relation
   *  \param the matrix to set
   *  \exception to be defined
   */
  inline void setC(SiconosMatrix &C)
  {
    this->C = C;
  };

  /** \fn void setD(SiconosMatrix)
   *  \brief allows to set the matrix D of the Relation
   *  \param the matrix to set
   *  \exception to be defined
   */
  inline void setD(SiconosMatrix &D)
  {
    this->D = D;
  };

  /** \fn void setE(SiconosMatrix)
   *  \brief allows to set the matrix E of the Relation
   *  \param the matrix to set
   *  \exception to be defined
   */
  inline void setE(SiconosMatrix &E)
  {
    this->E = E;
  };

  /** \fn void setA(SimpleVector&)
   *  \brief set vector A of the Relation
   *  \param SimpleVector& : new value of A
   *  \exception to be defined
   */
  inline void setA(/*SiconosVector*/SimpleVector &a)
  {
    this->a = a;
  };

  ////////////////////////////

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

  /** \fn void createRelation(LinearTIRXML * relationXML,
            SiconosMatrix* C, SiconosMatrix* D,
            SiconosMatrix* E, SiconosVector* a,
            Interaction * interaction)
   *  \brief allows to create the Relation with an xml file, or the needed data
   *  \param LinearTIRXML * : the XML object for this Relation
   *  \param SiconosMatrix* : the matrix C of this Relation
   *  \param SiconosMatrix* : the matrix D of this Relation
   *  \param SiconosMatrix* : the matrix E of this Relation
   *  \param SiconosVector* : the vector a of this Relation
   *  \exception RuntimeException
   */
  void createRelation(LinearTIRXML * relationXML,
                      SiconosMatrix* C = NULL, SiconosMatrix* D = NULL,
                      SiconosMatrix* E = NULL, SiconosVector* a = NULL);

  /** \fn LinearTIR* convert (Relation *r)
   *  \brief encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param Relation * : the relation which must be converted
   * \return a pointer on the relation if it is of the right type, NULL otherwise
   */
  static LinearTIR* convert(Relation *r);

protected:
  /** \fn void fillRelationWithRelationXML()
   *  \brief uses the RelationXML of the LinearTIR to fill the fields of this Relation
   */
  void fillRelationWithRelationXML();


private:
  /** a matrix specific to the LinearTIR \f$ y= C x + D \lambda + a \f$*/
  SiconosMatrix C;
  /** a matrix specific to the LinearTIR \f$ y= C x + D \lambda + a \f$*/
  SiconosMatrix D;
  /** a matrix specific to the LinearTIR \f$ r= B \lambda \f$*/
  SiconosMatrix E;
  /** a vector specific to the LinearTIR \f$ y= C x + D \lambda + a \f$*/
  /*SiconosVector*/
  SimpleVector a;

  //  /** the XML object linked to the LinearTIR to read XML data */
  //  LinearTIRXML *ltirelationxml;
};

#endif // LINEARTIRELATION_H
