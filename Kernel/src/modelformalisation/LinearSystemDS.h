#ifndef LINEARSYSTEMDS_H
#define LINEARSYSTEMDS_H

#include "LinearSystemDSXML.h"
#include "DynamicalSystem.h"

#include "SiconosMatrix.h"
//#include "SiconosVector.h"
#include "NewSiconosVector.h"
#include <iostream>
#include <vector>

using namespace std;


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

  /** \fn LinearSystemDS(DSXML*)
   *  \brief constructor with XML object of the LinearSystemDS
   *  \param DSXML* : the XML object corresponding
   */
  LinearSystemDS(DSXML*);
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

  /** \fn int getSizeU (void)
   *  \brief allow to get the size of the SiconosVector u
   *  \return the size of the SiconosVector u
   */
  //int getSizeU (void);

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
  /*SiconosVector*/
  SimpleVector* getUPtr(void);

  /** \fn SiconosVector* getFPtr (void)
   *  \brief get vector f
   *  \return SimpleVector* : pointer on f of the LinearSystemDS
   */
  /*SiconosVector*/
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
  inline void setU(/*SiconosVector*/SimpleVector &u)
  {
    this->u = u;
  };

  /** \fn void setF (SimpleVector&)
   *  \brief set vector f
   *  \param SimpleVector& : new value of f
   */
  inline void setF(/*SiconosVector*/SimpleVector &f)
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

  /** \fn void createDynamicalSystem(DSXML * nsdsXML, int number, int n,
                      SiconosVector* x0, NSDS * nsds)
   *  \brief allows to create the DynamicalSystem with an xml file, or the needed data
   *  \param DSXML * : the XML object for this DynamicalSystem
   *  \param int : the number for this DynamicalSystem
   *  \param int : the dimension of this DynamicalSystem
   *  \param SiconosVector* : the initial state of this DynamicalSystem
  //  *  \param NSDS * : The NSDS which contains this DynamicalSystem
   *  \exception RuntimeException
   */
  void createDynamicalSystem(DSXML * dsXML, int number = -1, int n = -1,
                             SiconosVector* x0 = NULL);//, NSDS * nsds = NULL);

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
  /*SiconosVector*/
  SimpleVector u;
  /** strength vector */
  /*SiconosVector*/
  SimpleVector f;

  /* contains the name of the plugin for u */
  string uFunctionName;
  /* contains the name of the plugin for f */
  string fFunctionName;

  /** class for manage plugin (open, close librairy...) */
  SiconosSharedLibrary cShared;

  //  /** \fn void (*vectorFieldPtr)(double time)
  //   * \brief compute the state
  //   * \param double : the time to make the computations
  //   */
  //  void (*vectorFieldPtr)(double time);

  //  /** \fn void (*computeJacobianPtr)(double time)
  //   * \brief compute the gradient of the state
  //   * \param double : the time to make the computations
  //   */
  //  void (*computeJacobianPtr)(double time);

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
//$Log: LinearSystemDS.h,v $
//Revision 1.39  2005/02/11 17:36:01  charlety
//
//_ little "inspection of code"
//_ basic getters and setters passed inline
//_ getters functions passed const
//
//Revision 1.38  2005/01/31 16:26:21  charlety
//
//_ Added a method named "convert" to classes which inherits from another. This is necessary for Python interface, in order to be able to use down-casting mechanism.
//
//Revision 1.37  2004/09/30 08:35:02  jbarbier
//- fonction of the formalisation : fill..With...XML and link... are now
//"protected" and no more "public"
//
//Revision 1.36  2004/09/22 11:16:28  charlety
//
//_ revision of Doxygen comments in modelformalisation
//
//Revision 1.35  2004/09/10 11:26:14  charlety
//
//_ Integration of the new version of the SiconosVector in the platform. the class simpleVector is used mostly to replace old SiconosVector. When a vector can be composite or simple, like the state of a dynamical system, a pointer on SiconosVector is used, and the vector is initialized simple or composite when the system is initialized.
//
//_ All the tests which worked with the previous version of the vector are OK with the new version.
//
//_ Example SICONOS and bouncingBall are OK
//
//_ some comments have still to be adapted to NewSiconosVector .
//
//_ This version of NewSiconosVector could be called 0.9. some details have to be fixed, it will be done before the end of September.
//
//Revision 1.34  2004/08/23 14:30:01  jbarbier
//- All the dynamical systems can be created in a comand program and added to a
//NSDS. The save is OK, but the creation of the boundary conditions is not yet
//finished.
//
//Revision 1.33  2004/08/20 15:26:44  jbarbier
//- creation of a Model and save in the XML is ok
//- creation of a NSDS and save in the XML is ok
//- creation of a NonLinearSystemDS and save in the XML is OK
//
//Revision 1.32  2004/08/13 11:26:58  jbarbier
//- function createNSDS complete
//
//- function createDynamicalSystem and createLinearSystemDS complete
//
//- function  createLagrangianNLDS in progress
//
//Revision 1.31  2004/08/12 11:55:15  jbarbier
//- new methods createModel, createNSDS, createStrategy, ...
//they now allow to make the link with upper objects of the platform
//it will be used for the creation of the platform without XML input file
//
//- the createModel method is finished but the attributes of the other objects
//of the platform are missing for the conctruction
//
//Revision 1.30  2004/07/29 14:25:37  jbarbier
