#ifndef LINEARDSIO_H
#define LINEARDSIO_H

#include "DSInputOutput.h"
#include "LinearDSIOXML.h"

/** \class LinearDSIO
 *  \brief Linear DSInputOutput
 *         { y = A.x
 *         { R = B.lambda
*  \author SICONOS Development Team - copyright INRIA
 *  \version 0.1
 *  \date Apr 27, 2004
 *
 *
 *
 */
class LinearDSIO : public DSInputOutput
{
public:

  /** \fn LinearDSIO();
   *  \brief Basic constructor
   */
  LinearDSIO();

  /** \fn LinearDSIO(LinearDSIOXML*)
   *  \brief constructor with XML object of the LinearDSIO
   *  \param LinearDSIOXML* : the XML object corresponding
   */
  LinearDSIO(DSInputOutputXML*);

  ~LinearDSIO();

  /** \fn void saveDSInputOutputToXML()
   *  \brief copy the data of the DSInputOutput to the XML tree
   */
  void saveDSInputOutputToXML();

  /** \fn void createDSInputOutput(DSInputOutputXML * dsioXML, int number = -1,
              SiconosMatrix *A = NULL, SiconosMatrix *B = NULL )
   *  \brief allows to create the DSInputOutput with an xml file, or the needed data
   *  \param DSInputOutputXML * : the XML object for this DSInputOutput
   *  \param int : the number of the DSInputOutput
   *  \param SiconosMatrix* : matrix A of the linear DSIO
   *  \param SiconosMatrix* : matrix B of the linear DSIO
   *  \exception RuntimeException
   */
  void createDSInputOutput(DSInputOutputXML * dsioXML, int number = -1,
                           SiconosMatrix *A = NULL, SiconosMatrix *B = NULL);


protected:
  /** \fn void fillDSInputOutputWithDSInputOutputXML()
   *  \brief uses the DSInputOutputXML of the LinearDSIO to fill the fields of this DSInputOutput
   */
  void fillDSInputOutputWithDSInputOutputXML();


private:
  SiconosMatrix A;
  SiconosMatrix B;
};

#endif // LINEARDSIO_H
