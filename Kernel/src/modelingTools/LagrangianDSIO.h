#ifndef LAGRANGIANDSIO_H
#define LAGRANGIANDSIO_H

#include "DSInputOutput.h"
#include "LagrangianDSIOXML.h"
#include "check.h"

extern std::string  DefaultComputeInput;
extern std::string  DefaultComputeOutput;


/** \class LagrangianDSIO
 *  \brief Lagrangian DSInputOutput
 *         { y = H(q, t)
 *         { R = G(lambda)
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

  //  /** \fn void createDSInputOutput(DSInputOutputXML * dsioXML)
  //   *  \brief allows to create the DSInputOutput with an xml file, or the needed data
  //   *  \param DSInputOutputXML * : the XML object for this DSInputOutput
  //   *  \exception RuntimeException
  //   */
  //  void createDSInputOutput(DSInputOutputXML * dsioXML, int number = -1,
  //                string computeInput=DefaultComputeInput,
  //                string computeOutput=DefaultComputeOutput);

private:

  /** class for manage plugin (open, close librairy...) */
  SiconosSharedLibrary cShared;

  //  /** \fn void (*computeJacobianPtr)(void);
  //   * \brief to be defined
  //   */
  //  void (*computeJacobianPtr)(int* sizeOfQ, double* qPtr, int* sizeOfY, double* jacobPtr);
  //
  //  void (*computeHPtr)(int* sizeOfQ, double* qPtr, int* sizeOfY, double* yPtr);

  /** \fn void (*computeOutputPtr)(double* xPtr, double* time, double* lambdaPtr, double* yPtr)
   *  \brief computes y
   *  \param double* xPtr : the pointer to the first element of the vector x
   *  \param double* time : the current time
   *  \param double* lambdaPtr : the pointer to the first element of the vector lambda
   *  \param double* yPtr : the pointer to the first element of the vector y (in-out parameter)
   */
  void (*computeOutputPtr)(double* xPtr, double* time, double* lambdaPtr, double* yPtr);

  /** \fn void (*computeInputPtr)(double* xPtr, double* time, double* lambdaPtr, double* rPtr)
   *  \brief computes r
   *  \param double* xPtr : the pointer to the first element of the vector x
   *  \param double* time : the current time
   *  \param double* lambdaPtr : the pointer to the first element of the vector lambda
   *  \param double* rPtr : the pointer to the first element of the vector r (in-out parameter)
   */
  void (*computeInputPtr)(double* xPtr, double* time, double* lambdaPtr, double* rPtr);
};

#endif
