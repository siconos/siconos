#ifndef LINEARSYSTEMDS_H
#define LINEARSYSTEMDS_H

#include "LinearSystemDSXML.h"
#include "DynamicalSystem.h"

#include "SiconosMatrix.h"
#include "NewSiconosVector.h"
#include <iostream>
#include <vector>

//using namespace std;


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

  /** \fn LinearSystemDS()
   *  \brief default constructor
   */
  LinearSystemDS();

  /** \fn LinearSystemDS(DSXML * nsdsXML)
   *  \brief allows to create the DynamicalSystem with an xml file, or the needed data
   *  \param DSXML * : the XML object for this DynamicalSystem
   *  \exception RuntimeException
   */
  LinearSystemDS(DSXML * dsXML);

  /** \fn LinearSystemDS(int number, int n,
                      SiconosVector* x0, NSDS * nsds)
   *  \brief allows to create the DynamicalSystem with an xml file, or the needed data
   *  \param int : the number for this DynamicalSystem
   *  \param int : the dimension of this DynamicalSystem
   *  \param SiconosVector* : the initial state of this DynamicalSystem
   *  \exception RuntimeException
   */
  LinearSystemDS(int number, int n, SiconosVector* x0);

  ~LinearSystemDS();

  // getter and setter

  /** \fn SiconosMatrix getA (void)
   *  \brief allow to get the size of the SiconosMatrix A
   *  \return the SiconosMatrix A of the LinearSystemDS
   */
  inline SiconosMatrix getA(void) const
  {
    return this->A;
  };

  /** \fn SiconosMatrix getB (void)
   *  \brief allow to get the size of the SiconosMatrix B
   *  \return the SiconosMatrix B of the LinearSystemDS
   */
  inline SiconosMatrix getB(void) const
  {
    return this->B;
  };

  /** \fn SiconosMatrix* getAPtr (void)
   *  \brief allow to get the size of the SiconosMatrix A
   *  \return the SiconosMatrix* A of the LinearSystemDS
   */
  SiconosMatrix* getAPtr(void);

  /** \fn SiconosMatrix* getBPtr (void)
   *  \brief allow to get the size of the SiconosMatrix B
   *  \return the SiconosMatrix* B of the LinearSystemDS
   */
  SiconosMatrix* getBPtr(void);

  /** \fn SimpleVector getU (void)
   *  \brief get the vector u
   *  \return SimpleVector : value of u
   */
  inline /*SiconosVector*/SimpleVector getU(void) const
  {
    return this->u;
  };

  /** \fn SimpleVector getF (void)
   *  \brief get vector f
   *  \return SimpleVector : value of f
   */
  inline /*SiconosVector*/SimpleVector getF(void) const
  {
    return this->f;
  };

  /** \fn SimpleVector* getUPtr (void)
   *  \brief allow to get the SiconosVector u
   *  \return the SiconosVector* u of the LinearSystemDS
   */
  SimpleVector* getUPtr(void);

  /** \fn SiconosVector* getFPtr (void)
   *  \brief get vector f
   *  \return SimpleVector* : pointer on f of the LinearSystemDS
   */
  SimpleVector* getFPtr(void);


  /** \fn void setA (SiconosMatrix)
   *  \brief allow to set the SiconosMatrix A
   *  \param a SiconosMatrix to set A
   */
  inline void setA(SiconosMatrix &A)
  {
    this->A = A;
  };

  /** \fn void setB (SiconosMatrix)
   *  \brief allow to set the SiconosMatrix B
   *  \param a SiconosMatrix to set B
   */
  inline void setB(SiconosMatrix &B)
  {
    this->B = B;
  };

  /** \fn void setU (SimpleVector&)
   *  \brief set vector u
   *  \param SimpleVector& : new value of u
   */
  inline void setU(SimpleVector &u)
  {
    this->u = u;
  };

  /** \fn void setF (SimpleVector&)
   *  \brief set vector f
   *  \param SimpleVector& : new value of f
   */
  inline void setF(SimpleVector &f)
  {
    this->f = f;
  };

  ////////////////////////////

  /** \fn void setComputeUFunction(string libPath, string functionName)
   *  \brief allow to set a specified function to compute the vector U
   *  \param string : the complete path to the plugin
   *  \param string : the function name to use in this plugin
   *  \exception SiconosSharedLibraryException
   */
  void setComputeUFunction(string pluginPath, string functionName);

  /** \fn void setComputeFFunction(string libPath, string functionName);
   *  \brief allow to set a specified function to compute the vector F
   *  \param string : the complete path to the plugin
   *  \param string : the function name to use in this plugin
   *  \exception SiconosSharedLibraryException
   */
  void setComputeFFunction(string pluginPath, string functionName);


  /** \fn void computeF(double time)
   *  \brief default function to compute vector F
   *  \exception RuntimeException
   */
  void computeF(double time);

  /** \fn void computeU(double time)
   *  \brief default function to compute vector U
   *  \exception RuntimeException
   */
  void computeU(double time);


  ////////////////////

  /** \fn void saveDSToXML()
   *  \brief copy the data of the DS to the XML tree
   *  \exception RuntimeException
   */
  void saveDSToXML();

  /** \fn void display()
   *  \brief print the data to the screen
   */
  void display() const;

  /** \fn LinearSystemDS* convert (DynamicalSystem* ds)
   *  \brief encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param DynamicalSystem* : the system which must be converted
   * \return a pointer on the dynamical system if it is of the right type, NULL otherwise
   */
  static LinearSystemDS* convert(DynamicalSystem* ds);

protected:
  /** \fn void fillDSWithDSXML()
   *  \brief overload of the function for a LinearSystemDS
   *  \exception RuntimeException
   */
  void fillDSWithDSXML();


private:
  /** matrix specific to the LinearSystemDS \f$ A \in R^{n \times n}  \f$*/
  SiconosMatrix A;
  /** matrix specific to the LinearSystemDS \f$ B \in R^{n \times uSize}  \f$ */
  SiconosMatrix B;
  /** size of vector u */
  int uSize;
  /** vector specific to the LinearSystemDS */
  SimpleVector u;
  /** strength vector */
  SimpleVector f;

  /* contains the name of the plugin for u */
  string uFunctionName;
  /* contains the name of the plugin for f */
  string fFunctionName;

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

  /** \fn void init()
   *  \brief initialise value of a LinearSystemDS
   */
  void init();
};

#endif // LINEARSYSTEMDS_H
