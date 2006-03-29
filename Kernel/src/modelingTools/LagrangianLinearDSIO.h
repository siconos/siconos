/* Siconos-Kernel version 1.1.4, Copyright INRIA 2005-2006.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
*/
#ifndef LAGRANGIANLINEARDSIO_H
#define LAGRANGIANLINEARDSIO_H

#include "LagrangianDSIO.h"
#include "LagrangianLinearDSIOXML.h"

/** \class LagrangianLinearDSIO
 *  \brief Lagrangian DSInputOutput
 *         { y = H.q + b
 *         { R = Ht
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 1.1.4.
 *  \date 17/01/2005
 *
 *
 */
class LagrangianLinearDSIO : public LagrangianDSIO
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


  /** \fn SimpleMatrix getH(void)
   *  \brief getter of the SimpleMatrix H
   *  \return a pointer on the SimpleMatrix H
   */
  inline SimpleMatrix getH(void) const
  {
    return *H;
  } ;

  /** \fn SimpleVector getB(void)
   *  \brief getter of the SiconosVector b
   *  \return SimpleVector : value of b
   */
  inline SimpleVector getB(void) const
  {
    return *b;
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
  inline void setH(const SiconosMatrix &newH)
  {
    *H = newH;
  };

  /** \fn void setH(SimpleVector&)
   *  \brief set the vector b
   *  \param SimpleVector& : new value of b
   */
  inline void setB(SimpleVector& newB)
  {
    *b = newB;
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
  SiconosMatrix* H;

  /** a specific vector to the LagrangianLinearDSIO */
  SimpleVector* b;
};

#endif // LAGRANGIANLINEARDSIO_H
