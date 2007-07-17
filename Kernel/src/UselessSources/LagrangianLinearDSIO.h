/* Siconos-Kernel version 2.1.1, Copyright INRIA 2005-2006.
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
/*! \file LagrangianLinearDSIO.h

*/
#ifndef LAGRANGIANLINEARDSIO_H
#define LAGRANGIANLINEARDSIO_H

#include "LagrangianDSIO.h"
#include "LagrangianLinearDSIOXML.h"

//! Lagrangian DSInputOutput
/**         { y = H.q + b
 *         { R = Ht
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 2.1.1.
 *  \date 17/01/2005
 *
 *
 */
class LagrangianLinearDSIO : public LagrangianDSIO
{
public:
  /** Default constructor
  */
  LagrangianLinearDSIO();

  /** constructor with XML object of the parent class DSInputOutput
  *  \param DSInputOutputXML* : the XML object corresponding
  */
  LagrangianLinearDSIO(DSInputOutputXML*);

  virtual ~LagrangianLinearDSIO();


  /** getter of the SimpleMatrix H
  *  \return a pointer on the SimpleMatrix H
  */
  inline SimpleMatrix getH(void) const
  {
    return *H;
  } ;

  /** getter of the SiconosVector b
  *  \return SimpleVector : value of b
  */
  inline SimpleVector getB(void) const
  {
    return *b;
  };

  /** getter of the SiconosMatrix* h
  *  \return a pointer on the SiconosMatrix* h
  */
  SiconosMatrix* getHPtr(void);

  /** getter of the SiconosVector* b
  *  \return a pointer on the SiconosVector b
  */
  SiconosVector* getBPtr(void);

  /** setter on the SiconosMatrix h
  *  \param a SiconosMatrix to set h
  */
  inline void setH(const SiconosMatrix &newH)
  {
    *H = newH;
  };

  /** set the vector b
  *  \param SimpleVector& : new value of b
  */
  inline void setB(SimpleVector& newB)
  {
    *b = newB;
  };

  /** default function to compute y for the free state
  *  \param double : current time
  *  \exception RuntimeException
  */
  void computeFreeOutput(double time);

  /** default function to compute y
  *  \param double : current time
  *  \exception RuntimeException
  */
  void computeOutput(double time);

  /** default function to compute lambda
  *  \param double : current time
  *  \exception RuntimeException
  */
  void computeInput(double time);

  /** copy the data of the DSInputOutput to the XML tree
  *  \exception RuntimeException
  */
  void saveDSInputOutputToXML();

  /** print the data to the screen
  */
  void display() const;

  /** allows to create the DSInputOutput with an xml file, or the needed data
   *  \param DSInputOutputXML * : the XML object for this DSInputOutput
   *  \param SiconosMatrix* : the matrix H of this DSInputOutput
   *  \param SiconosVector* : the vector b of this DSInputOutput
   *  \exception RuntimeException
   */
  void createDSInputOutput(DSInputOutputXML * dsioXML,
                           SiconosMatrix* H = NULL, SiconosVector* b = NULL);

  /** encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param Relation * : the DSInputOutput which must be converted
   * \return a pointer on the DSInputOutput if it is of the right type, NULL otherwise
   */
  static LagrangianLinearDSIO* convert(DSInputOutput *r);

protected:
  /** uses the DSInputOutputXML of the LagrangianLinearDSIO to fill the fields of this DSInputOutput
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
