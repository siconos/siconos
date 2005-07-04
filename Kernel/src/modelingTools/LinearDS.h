#ifndef LINEARDS_H
#define LINEARDS_H

#include "LinearDSXML.h"
#include "DynamicalSystem.h"
#include "SiconosMatrix.h"
#include "NewSiconosVector.h"
#include "check.h"
#include <iostream>
#include <vector>

class LinearDSXML;

/** \class LinearDS
 *  \brief First order linear systems - Inherits from DynamicalSystems
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 1.0
 *  \date (Creation) Apr 29, 2004
 *
 *
 *  This class represents first order linear systems of the form:
 * \f[
 * \dot x = A(t)x(t)+B(t)u(t)+f(t)+r,
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
 *
 **/

class LinearDS : public DynamicalSystem
{
public:

  /** \fn LinearDS(DSXML * nsdsXML)
   *  \brief xml constructor
   *  \param DSXML * : the XML object for this DynamicalSystem
   *  \exception RuntimeException
   */
  LinearDS(DSXML * dsXML);

  /** \fn LinearDS(int number, int n, SiconosVector* x0, NSDS * nsds)
   *  \brief constructor from a set of data
   *  \param int : reference number of this DynamicalSystem
   *  \param int : dimension of this DynamicalSystem
   *  \param SiconosVector : the initial state of this DynamicalSystem
   *  \param string: plugin path for A (optional)
   *  \param string: plugin function name for A (optional)
   *  \exception RuntimeException
   */
  LinearDS(const int&, const unsigned int&, const SiconosVector&,
           const std::string& = "BasicPlugin.so", const std::string& = "computeA");

  /** \fn LinearDS()
   *  \brief constructor from a set of data
   *  \param int : reference number of the DynamicalSystem
   *  \param SiconosVector : the initial state of this DynamicalSystem
   *  \param SiconosMatrix : matrix A
   *  \exception RuntimeException
   */
  LinearDS(const int& newNumber, const SiconosVector& newX0,
           const SiconosMatrix& newA);

  /** \fn ~LinearDS()
   *  \brief destructor */
  ~LinearDS();

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
    isPlugin[0] = false;
  }

  /** \fn void setAPtr(SiconosMatrix* newPtr)
   *  \brief set A to pointer newPtr
   *  \param SiconosMatrix * newPtr
   */
  void setAPtr(SiconosMatrix *);

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
  void setF(const SimpleVector&);

  /** \fn void setFPtr(SimpleVector* newPtr)
   *  \brief set F to pointer newPtr
   *  \param SimpleVector * newPtr
   */
  void setFPtr(SimpleVector *);

  // uSize

  /** \fn const int getUSize() const
   *  \brief get the value of uSize
   *  \return the value of uSize
   */
  inline const unsigned int getUSize() const
  {
    return uSize;
  };

  /** \fn void setUSize(const int&)
   *  \brief set uSize AND allocate memory for u and B
   *  \param int uSize : the value to set uSize
   */
  void setUSize(const unsigned int&);

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
  void setU(const SimpleVector&);

  /** \fn void setUPtr(SimpleVector* newPtr)
   *  \brief set U to pointer newPtr
   *  \param SimpleVector * newPtr
   */
  void setUPtr(SimpleVector *);

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
  void setB(const SiconosMatrix&);

  /** \fn void setBPtr(SiconosMatrix* newPtr)
   *  \brief set B to pointer newPtr
   *  \param SiconosMatrix * newPtr
   */
  void setBPtr(SiconosMatrix *);

  // --- plugins related functions

  /** \fn void setComputeAFunction(const string& libPath,const string& functionName)
   *  \brief set a specified function to compute the matrix A
   *  \param string : the complete path to the plugin
   *  \param string : the function name to use in this plugin
   *  \exception SiconosSharedLibraryException
   */
  void setComputeAFunction(const std::string &, const std::string &);

  /** \fn void setComputeFFunction(const string& libPath,const string& functionName);
   *  \brief set a specified function to compute the vector F
   *  \param string : the complete path to the plugin
   *  \param string : the function name to use in this plugin
   *  \exception SiconosSharedLibraryException
   */
  void setComputeFFunction(const std::string &, const std::string &);

  /** \fn void setComputeUFunction(const string& libPath,const string& functionName)
   *  \brief set a specified function to compute the vector U
   *  \param string : the complete path to the plugin
   *  \param string : the function name to use in this plugin
   *  \exception SiconosSharedLibraryException
   */
  void setComputeUFunction(const std::string &, const std::string &);

  /** \fn void setComputeBFunction(const string& libPath,const string& functionName);
   *  \brief set a specified function to compute the matrix B
   *  \param string : the complete path to the plugin
   *  \param string : the function name to use in this plugin
   *  \exception SiconosSharedLibraryException
   */
  void setComputeBFunction(const std::string &, const std::string &);

  /** \fn void computeA(const double& time)
   *  \brief default function to compute matrix A
   *  \exception RuntimeException
   */
  void computeA(const double&);

  /** \fn void computeF(const double& time)
   *  \brief default function to compute vector F
   *  \exception RuntimeException
   */
  void computeF(const double&);

  /** \fn void computeU(const double& time)
   *  \brief default function to compute vector U
   *  \exception RuntimeException
   */
  void computeU(const double&);

  /** \fn void computeB(const double& time)
   *  \brief default function to compute matrix B
   *  \exception RuntimeException
   */
  void computeB(const double&);

  // --- xml related functions ---

  /** \fn void saveDSToXML()
   *  \brief copy the data of the DS into the XML tree
   *  \exception RuntimeException
   */
  void saveDSToXML();

  /** \fn void display()
   *  \brief data display on screen
   */
  void display() const;

  /** \fn LinearDS* convert (DynamicalSystem* ds)
   *  \brief encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param DynamicalSystem* : the system which must be converted
   * \return a pointer on the dynamical system if it is of the right type, NULL otherwise
   */
  static LinearDS* convert(DynamicalSystem* ds);

private:

  /** \fn LinearDS()
   *  \brief default constructor
   */
  LinearDS();

  /** matrix specific to the LinearDS \f$ A \in R^{n \times n}  \f$*/
  SiconosMatrix *A;
  /** strength vector */
  SimpleVector *f;
  /** size of vector u = number of columns in B*/
  unsigned int uSize;
  /** vector specific to the LinearDS */
  SimpleVector *u;
  /** matrix specific to the LinearDS \f$ B \in R^{n \times uSize}  \f$ */
  SiconosMatrix *B;

  /* contains the name of the plugin for A */
  std::string  AFunctionName;
  /* contains the name of the plugin for f */
  std::string  fFunctionName;
  /* contains the name of the plugin for u */
  std::string  uFunctionName;
  /* contains the name of the plugin for B */
  std::string  BFunctionName;

  /** class for manage plugin (open, close librairy...) */
  SiconosSharedLibrary cShared;

  /** \fn void (*computeAPtr)(int sizeOfA, double* fPtr,double time)
   *  \brief compute matrix A
   *  \param double : the time to make the computations
   */
  void (*computeAPtr)(unsigned int* sizeOfA, double* APtr, const double* time);

  /** \fn void (*computeFPtr)(int sizeOfF, double* fPtr,double time)
   *  \brief compute vector F
   *  \param double : the time to make the computations
   */
  void (*computeFPtr)(unsigned int* sizeOfF, double* fPtr, const double* time);

  /** \fn void (*computeUPtr)(int sizeOfU, double* uPtr,double time)
   *  \brief compute vector U
   *  \param double : the time to make the computations
   */
  void (*computeUPtr)(unsigned int* sizeOfU, double* uPtr, const double* time);

  /** \fn void (*computeUPtr)(int rowsOfB, int colOfB, double* BPtr,double time)
   *  \brief compute matrix B
   *  \param double : the time to make the computations
   */
  void (*computeBPtr)(unsigned int* rowsOfB, unsigned int* colOfB, double* BPtr, const double* time);

  /** vector of bool to check if A, f, u, B (in this order!) are loaded from a plugin or not */
  std::vector<bool> isPlugin;

  /** Flags to know if pointers have been allocated inside constructors or not */

  bool isAAllocatedIn;
  bool isFAllocatedIn;
  bool isUAllocatedIn;
  bool isBAllocatedIn;
};

#endif // LINEARDS_H
