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
 * \dot x = A(t)x(t)+E(t)u(t)+b(t)+r,
 *  x(t_0)=x_0
 * \f]
 * where
 *    - \f$x \in R^{n} \f$ is the state,
 *    - \f$r \in R^{n} \f$  the input due to the Non Smooth Interaction.
 *
 *  The  VectorField is specialized by
 *    - \f$A \in R^{n\times n} \f$
 *    - \f$b \in R^{n} \f$
 *    - \f$u \in R^{uSize} \f$
 *    - \f$E \in R^{n\times uSize} \f$
 *
 * The "minimal" form is
 * \f[
 * \dot x = A(t)x(t),
 *  x(t_0)=x_0
 * \f]
 * and so A should always be specified.
 *
 *  \todo Automatically, specify the function of DynamicalSystem such as
 *          VectorField.
 *
 **/

class LinearDS : public DynamicalSystem
{
public:

  /** === CONSTRUCTORS/DESTRUCTOR === */

  /** \fn LinearDS(DSXML * nsdsXML)
   *  \brief xml constructor
   *  \param DSXML * : the XML object for this DynamicalSystem
   *  \param NonSmoothDynamicalSystem* (optional): the NSDS that owns this ds
   *  \exception RuntimeException
   */
  LinearDS(DSXML * dsXML, NonSmoothDynamicalSystem* = NULL);

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

  /** \fn LinearDS( const int& newNumber, const SiconosVector& newX0,
   *                const SiconosMatrix& newA)
   *  \brief constructor from a set of data
   *  \param int : reference number of the DynamicalSystem
   *  \param SiconosVector : the initial state of this DynamicalSystem
   *  \param SiconosMatrix : matrix A
   *  \exception RuntimeException
   */
  LinearDS(const int& newNumber, const SiconosVector& newX0,
           const SiconosMatrix& newA);

  /** \fn LinearDS(const DynamicalSystem &)
   *  \brief copy constructor
   *  \param a Dynamical system to copy
   */
  LinearDS(const DynamicalSystem &);

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

  // --- b ---

  /** \fn  const SimpleVector getB() const
   *  \brief get the value of b
   *  \return SimpleVector
   */
  inline const SimpleVector getB() const
  {
    return *b;
  }

  /** \fn SimpleVector* getBPtr() const
   *  \brief get b
   *  \return pointer on a SimpleVector
   */
  inline SimpleVector* getBPtr() const
  {
    return b;
  }

  /** \fn void setB (const SimpleVector& newValue)
   *  \brief set the value of b to newValue
   *  \param SimpleVector newValue
   */
  void setB(const SimpleVector&);

  /** \fn void setBPtr(SimpleVector* newPtr)
   *  \brief set b to pointer newPtr
   *  \param SimpleVector * newPtr
   */
  void setBPtr(SimpleVector *);

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

  // --- E ---
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
  void setE(const SiconosMatrix&);

  /** \fn void setEPtr(SiconosMatrix* newPtr)
   *  \brief set E to pointer newPtr
   *  \param SiconosMatrix * newPtr
   */
  void setEPtr(SiconosMatrix *);

  // --- plugins related functions

  /** \fn  std::vector<bool> getIsPlugin() const
   *  \brief get boolean vector that checks if members are loaded from plugin or not
   *  \return a vector of bool
   */
  inline const std::vector<bool> getIsPlugin() const
  {
    return isPlugin;
  }

  /** \fn  std::string getAFunctionName() const
  *  \brief get name of function that computes A (if A from plugin)
  *  \return a string
  */
  inline const std::string getAFunctionName() const
  {
    return AFunctionName;
  }

  /** \fn  std::string getBFunctionName() const
  *  \brief get name of function that computes b (if b from plugin)
  *  \return a string
  */
  inline const std::string getBFunctionName() const
  {
    return bFunctionName;
  }

  /** \fn  std::string getUFunctionName() const
  *  \brief get name of function that computes u (if u from plugin)
  *  \return a string
  */
  inline const std::string getUFunctionName() const
  {
    return uFunctionName;
  }

  /** \fn  std::string getEFunctionName() const
  *  \brief get name of function that computes E (if E from plugin)
  *  \return a string
  */
  inline const std::string getEFunctionName() const
  {
    return EFunctionName;
  }

  /** \fn void setComputeAFunction(const string& libPath,const string& functionName)
   *  \brief set a specified function to compute the matrix A
   *  \param string : the complete path to the plugin
   *  \param string : the function name to use in this plugin
   *  \exception SiconosSharedLibraryException
   */
  void setComputeAFunction(const std::string &, const std::string &);

  /** \fn void setComputeBFunction(const string& libPath,const string& functionName);
   *  \brief set a specified function to compute the vector b
   *  \param string : the complete path to the plugin
   *  \param string : the function name to use in this plugin
   *  \exception SiconosSharedLibraryException
   */
  void setComputeBFunction(const std::string &, const std::string &);

  /** \fn void setComputeUFunction(const string& libPath,const string& functionName)
   *  \brief set a specified function to compute the vector u
   *  \param string : the complete path to the plugin
   *  \param string : the function name to use in this plugin
   *  \exception SiconosSharedLibraryException
   */
  void setComputeUFunction(const std::string &, const std::string &);

  /** \fn void setComputeEFunction(const string& libPath,const string& functionName);
   *  \brief set a specified function to compute the matrix E
   *  \param string : the complete path to the plugin
   *  \param string : the function name to use in this plugin
   *  \exception SiconosSharedLibraryException
   */
  void setComputeEFunction(const std::string &, const std::string &);

  /** \fn void computeA(const double& time)
   *  \brief default function to compute matrix A
   *  \exception RuntimeException
   */
  void computeA(const double&);

  /** \fn void computeB(const double& time)
   *  \brief default function to compute vector b
   *  \exception RuntimeException
   */
  void computeB(const double&);

  /** \fn void computeU(const double& time)
   *  \brief default function to compute vector U
   *  \exception RuntimeException
   */
  void computeU(const double&);

  /** \fn void computeE(const double& time)
   *  \brief default function to compute matrix E
   *  \exception RuntimeException
   */
  void computeE(const double&);

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
  SimpleVector *b;
  /** size of vector u = number of columns in B*/
  unsigned int uSize;
  /** vector specific to the LinearDS */
  SimpleVector *u;
  /** matrix specific to the LinearDS \f$ B \in R^{n \times uSize}  \f$ */
  SiconosMatrix *E;

  /* contains the name of the plugin for A */
  std::string  AFunctionName;
  /* contains the name of the plugin for b */
  std::string  bFunctionName;
  /* contains the name of the plugin for u */
  std::string  uFunctionName;
  /* contains the name of the plugin for E */
  std::string  EFunctionName;

  /** class for plugin management (open, close library...) */
  SiconosSharedLibrary cShared;

  /** \fn void (*computeAPtr)(int sizeOfA, double* APtr,double time)
   *  \brief compute matrix A
   *  \param double : the time to make the computations
   */
  void (*computeAPtr)(unsigned int* sizeOfA, double* APtr, const double* time);

  /** \fn void (*computeBPtr)(int sizeOfB, double* bPtr,double time)
   *  \brief compute vector b
   *  \param double : the time to make the computations
   */
  void (*computeBPtr)(unsigned int* sizeOfB, double* bPtr, const double* time);

  /** \fn void (*computeUPtr)(int sizeOfU, double* uPtr,double time)
   *  \brief compute vector U
   *  \param double : the time to make the computations
   */
  void (*computeUPtr)(unsigned int* sizeOfU, double* uPtr, const double* time);

  /** \fn void (*computeEPtr)(int rowsOfE, int colOfE, double* EPtr,double time)
   *  \brief compute matrix E
   *  \param double : the time to make the computations
   */
  void (*computeEPtr)(unsigned int* rowsOfE, unsigned int* colOfE, double* EPtr, const double* time);

  /** vector of bool to check if A, b, u, E (in this order!) are loaded from a plugin or not */
  std::vector<bool> isPlugin;

  /** Flags to know if pointers have been allocated inside constructors or not */

  bool isAAllocatedIn;
  bool isBAllocatedIn;
  bool isUAllocatedIn;
  bool isEAllocatedIn;
};

#endif // LINEARDS_H
