#ifndef LAGRANGIANDSIO_H
#define LAGRANGIANDSIO_H

#include "DSInputOutput.h"
#include "LagrangianDSIOXML.h"

/** \class LagrangianDSIO
 *  \brief Lagrangian DSInputOutput
*  \author SICONOS Development Team - copyright INRIA
 *  \version 0.1
 *  \date 17/01/2005
 *
 *
 */
class LagrangianDSIO : public DSInputOutput
{
public:

  /** \fn LagrangianDSIO(void);
   * \brief default constructor
   */
  LagrangianDSIO();
  ~LagrangianDSIO();

  //  /** \fn void saveDSInputOutputToXML()
  //   *  \brief copy the data of the DSInputOutput to the XML tree
  //   */
  //  void saveDSInputOutputToXML();

  /** \fn void createDSInputOutput(DSInputOutputXML * dsioXML)
   *  \brief allows to create the DSInputOutput with an xml file, or the needed data
   *  \param DSInputOutputXML * : the XML object for this DSInputOutput
   *  \exception RuntimeException
   */
  void createDSInputOutput(DSInputOutputXML * dsioXML, int number = -1,
                           SiconosMatrix *H = NULL);

private:

  //  /** class for manage plugin (open, close librairy...) */
  //  SiconosSharedLibrary cShared;
  //
  //  /** \fn void (*computeJacobianPtr)(void);
  //   * \brief to be defined
  //   */
  //
  //  void (*computeJacobianPtr)(int* sizeOfQ, double* qPtr, int* sizeOfY, double* jacobPtr);
  //
  //  void (*computeHPtr)(int* sizeOfQ, double* qPtr, int* sizeOfY, double* yPtr);

};

#endif
