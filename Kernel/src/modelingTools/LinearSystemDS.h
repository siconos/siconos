#ifndef LINEARSYSTEMDS_H
#define LINEARSYSTEMDS_H

#include "LinearSystemDSXML.h"
#include "DynamicalSystem.h"
#include "SiconosMatrix.h"
#include "NewSiconosVector.h"
#include "check.h"
#include <iostream>
#include <vector>

class LinearSystemDSXML;

/** \class LinearSystemDS
 *  \brief class of dynamic systems, inherited of DynamicalSystem
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 1.0
 *  \date (Creation) Apr 29, 2004
 *
 *
 *  The class DynamicalSystem allows to define and compute a generic n-dimensional
 * linear dynamical system of the form :
 * \f[
 * \dot x = Ax+Bu+f+r,
 *  x(t_0)=x_0
 * \f]
 * where
 *    - \f$x \in R^{n} \f$ is the state,
 *    - \f$r \in R^{n} \f$  the input due to the Non Smooth Interaction.
 *
 *  The  VectorField is specialized by
 *    - \f$A \in R^{n\times n} \f$
 *    - \f$f \in R^{n} \f$
 *    - \f$u \in R^{uSize} \f$
 *    - \f$B \in R^{n\times uSize} \f$
 *
 *  \todo Automatically, specify the function of DynamicalSystem such as
 *          VectorField.
 **/

class LinearSystemDS : public DynamicalSystem
{
public:

  /** \fn LinearSystemDS(DSXML * nsdsXML)
   *  \brief create the DynamicalSystem with an xml file
   *  \param DSXML * : the XML object for this DynamicalSystem
   *  \exception RuntimeException
   */
  LinearSystemDS(DSXML * dsXML);

  /** \fn LinearSystemDS(int number, int n, SiconosVector* x0, NSDS * nsds)
   *  \brief create the DynamicalSystem with a minimum set of data
   *  \param int : the number for this DynamicalSystem
   *  \param int : the dimension of this DynamicalSystem
   *  \param SiconosVector* : the initial state of this DynamicalSystem
   *  \exception RuntimeException
   */
  LinearSystemDS(int number, int n, SiconosVector* x0);

  ~LinearSystemDS();

  // --- getter and setter ---

  // --- A ---
  /** \fn  const SiconosMatrix getA() const
   *  \brief get the value of A
   *  \return SiconosMatrix
   */
  inline const SiconosMatrix getA() const
  {
    return *A;
  }

  /** \fn SiconosMatrix* getAPtr() const
   *  \brief get A
   *  \return pointer on a SiconosMatrix
   */
  inline SiconosMatrix* getAPtr() const
  {
    return A;
  }

  /** \fn void setA (const SiconosMatrix& newValue)
   *  \brief set the value of A to newValue
   *  \param SiconosMatrix newValue
   */
  inline void setA(const SiconosMatrix& newValue)
  {
    *A = newValue;
  }

  /** \fn void setAPtr(SiconosMatrix* newPtr)
   *  \brief set A to pointer newPtr
   *  \param SiconosMatrix * newPtr
   */
  inline void setAPtr(SiconosMatrix *newPtr)
  {
    delete A ;
    A = newPtr;
  }

  // --- B ---
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
  inline void setB(const SiconosMatrix& newValue)
  {
    *B = newValue;
  }

  /** \fn void setBPtr(SiconosMatrix* newPtr)
   *  \brief set B to pointer newPtr
   *  \param SiconosMatrix * newPtr
   */
  inline void setBPtr(SiconosMatrix *newPtr)
  {
    delete B;
    B = newPtr;
  }

  // --- u ---

  /** \fn  const SimpleVector getU() const
  *  \brief get the value of u
  *  \return SimpleVector
  */
  inline const SimpleVector getU() const
  {
    return *u;
  }

  /** \fn SimpleVector* getUPtr() const
   *  \brief get u
   *  \return pointer on a SimpleVector
   */
  inline SimpleVector* getUPtr() const
  {
    return u;
  }

  /** \fn void setU (const SimpleVector& newValue)
   *  \brief set the value of u to newValue
   *  \param SimpleVector newValue
   */
  inline void setU(const SimpleVector& newValue)
  {
    *u = newValue;
  }

  /** \fn void setUPtr(SimpleVector* newPtr)
   *  \brief set U to pointer newPtr
   *  \param SimpleVector * newPtr
   */
  inline void setUPtr(SimpleVector *newPtr)
  {
    delete u;
    u = newPtr;
  }

  // --- f ---

  /** \fn  const SimpleVector getF() const
   *  \brief get the value of f
   *  \return SimpleVector
   */
  inline const SimpleVector getF() const
  {
    return *f;
  }

  /** \fn SimpleVector* getFPtr() const
   *  \brief get f
   *  \return pointer on a SimpleVector
   */
  inline SimpleVector* getFPtr() const
  {
    return f;
  }

  /** \fn void setF (const SimpleVector& newValue)
   *  \brief set the value of f to newValue
   *  \param SimpleVector newValue
   */
  inline void setF(const SimpleVector& newValue)
  {
    *f = newValue;
  }

  /** \fn void setFPtr(SimpleVector* newPtr)
   *  \brief set F to pointer newPtr
   *  \param SimpleVector * newPtr
   */
  inline void setFPtr(SimpleVector *newPtr)
  {
    delete f;
    f = newPtr;
  }

  // --- plugins related functions

  /** \fn void setComputeUFunction(const string& libPath,const string& functionName)
   *  \brief set a specified function to compute the vector U
   *  \param string : the complete path to the plugin
   *  \param string : the function name to use in this plugin
   *  \exception SiconosSharedLibraryException
   */
  void setComputeUFunction(const std::string & pluginPath, const std::string & functionName);

  /** \fn void setComputeFFunction(const string& libPath,const string& functionName);
   *  \brief set a specified function to compute the vector F
   *  \param string : the complete path to the plugin
   *  \param string : the function name to use in this plugin
   *  \exception SiconosSharedLibraryException
   */
  void setComputeFFunction(const std::string & pluginPath, const std::string & functionName);

  /** \fn void computeF(const double& time)
   *  \brief default function to compute vector F
   *  \exception RuntimeException
   */
  void computeF(const double& time);

  /** \fn void computeU(const double& time)
   *  \brief default function to compute vector U
   *  \exception RuntimeException
   */
  void computeU(const double& time);

  ////////////////////

  /** \fn void saveDSToXML()
   *  \brief copy the data of the DS into the XML tree
   *  \exception RuntimeException
   */
  void saveDSToXML();

  /** \fn void display()
   *  \brief data display on screen
   */
  void display() const;

  /** \fn LinearSystemDS* convert (DynamicalSystem* ds)
   *  \brief encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param DynamicalSystem* : the system which must be converted
   * \return a pointer on the dynamical system if it is of the right type, NULL otherwise
   */
  static LinearSystemDS* convert(DynamicalSystem* ds);

private:

  /** \fn LinearSystemDS()
   *  \brief default constructor
   */
  LinearSystemDS();

  /** matrix specific to the LinearSystemDS \f$ A \in R^{n \times n}  \f$*/
  SiconosMatrix *A;
  /** matrix specific to the LinearSystemDS \f$ B \in R^{n \times uSize}  \f$ */
  SiconosMatrix *B;
  /** size of vector u */
  int uSize;
  /** vector specific to the LinearSystemDS */
  SimpleVector *u;
  /** strength vector */
  SimpleVector *f;

  /* contains the name of the plugin for u */
  std::string  uFunctionName;
  /* contains the name of the plugin for f */
  std::string  fFunctionName;

  /** class for manage plugin (open, close librairy...) */
  SiconosSharedLibrary cShared;

  /** \fn void (*computeFPtr)(int sizeOfF, double* fPtr,double time)
   *  \brief compute vector F
   *  \param double : the time to make the computations
   */
  void (*computeFPtr)(int* sizeOfF, double* fPtr, double* time);

  /** \fn void (*computeUPtr)(int sizeOfU, double* uPtr,double time)
   *  \brief compute vector U
   *  \param double : the time to make the computations
   */
  void (*computeUPtr)(int* sizeOfU, double* uPtr, double* time);

private :

};

#endif // LINEARSYSTEMDS_H
