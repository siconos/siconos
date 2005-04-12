#ifndef LAGRANGIANLINEARDSIO_H
#define LAGRANGIANLINEARDSIO_H

#include "DSInputOutput.h"
#include "LagrangianLinearDSIOXML.h"
#include "SimpleVector.h"
#include "CompositeVector.h"

/** \class LagrangianLinearDSIO
 *  \brief Lagrangian DSInputOutput
 *         { y = H.q + b
 *         { R = Ht
*  \author SICONOS Development Team - copyright INRIA
 *  \version 0.1
 *  \date 17/01/2005
 *
 *
 */
class LagrangianLinearDSIO : public DSInputOutput
{
public:
  /** \fn LagrangianLinearDSIO()
   *  \brief Default constructor
   */
  LagrangianLinearDSIO();

  /** \fn LagrangianLinearDSIO(DSInputOutputXML*)
   *  \brief constructor with XML object of the parent class DSInputOutput
   *  \param DSInputOutputXML* : the XML object corresponding
   */
  LagrangianLinearDSIO(DSInputOutputXML*);

  virtual ~LagrangianLinearDSIO();


  /** \fn SiconosMatrix getH(void)
   *  \brief getter of the SiconosMatrix H
   *  \return a pointer on the SiconosMatrix H
   */
  inline SiconosMatrix getH(void) const
  {
    return this->H;
  } ;

  /** \fn SimpleVector getB(void)
   *  \brief getter of the SiconosVector b
   *  \return SimpleVector : value of b
   */
  inline SimpleVector getB(void) const
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
    this->H = H;
  };

  /** \fn void setH(SimpleVector&)
   *  \brief set the vector b
   *  \param SimpleVector& : new value of b
   */
  inline void setB(SimpleVector& B)
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

  /** \fn void saveDSInputOutputToXML()
   *  \brief copy the data of the DSInputOutput to the XML tree
   *  \exception RuntimeException
   */
  void saveDSInputOutputToXML();

  /** \fn void display()
   *  \brief print the data to the screen
   */
  void display() const;

  /** \fn void createDSInputOutput(DSInputOutputXML * dsioXML,
            SiconosMatrix* H, SiconosVector* b)
   *  \brief allows to create the DSInputOutput with an xml file, or the needed data
   *  \param DSInputOutputXML * : the XML object for this DSInputOutput
   *  \param SiconosMatrix* : the matrix H of this DSInputOutput
   *  \param SiconosVector* : the vector b of this DSInputOutput
   *  \exception RuntimeException
   */
  void createDSInputOutput(DSInputOutputXML * dsioXML,
                           SiconosMatrix* H = NULL, SiconosVector* b = NULL);

  /** \fn LagrangianLinearDSIO* convert (DSInputOutput *r)
   *  \brief encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param Relation * : the DSInputOutput which must be converted
   * \return a pointer on the DSInputOutput if it is of the right type, NULL otherwise
   */
  static LagrangianLinearDSIO* convert(DSInputOutput *r);

protected:
  /** \fn void fillDSInputOutputWithDSInputOutputXML()
   *  \brief uses the DSInputOutputXML of the LagrangianLinearDSIO to fill the fields of this DSInputOutput
   *  \exception RuntimeException
   */
  void fillDSInputOutputWithDSInputOutputXML();

private:
  /** a specific matrix to the LagrangianLinearDSIO */
  SiconosMatrix H;

  /** a specific vector to the LagrangianLinearDSIO */
  SimpleVector b;
};

#endif // LAGRANGIANLINEARDSIO_H
