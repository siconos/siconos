//$Id: LagrangianTIDS.h,v 1.32 2005/02/11 17:36:00 charlety Exp $
#ifndef LAGRANGIANTIDS_H
#define LAGRANGIANTIDS_H

#include "LagrangianNLDS.h"
#include "LagrangianTIDSXML.h"

#include "SiconosMatrix.h"
//#include "SiconosVector.h"
#include "NewSiconosVector.h"

#include <iostream>
#include <vector>

using namespace std;


class LagrangianTIDSXML;

/** \class LagrangianTIDS
 *  \brief class of Lagrangian invariant time systems, inherited of LagrangianNLDS
 *  \author V. ACARY, JB. CHARLETY
 *  \version 1.0
 *  \date (Creation) Apr 29, 2004
 *
 * $Date: 2005/02/11 17:36:00 $
 * $Revision: 1.32 $
 * $Author: charlety $
 * $Source: /CVS/Siconos/SICONOS/src/modelformalisation/LagrangianTIDS.h,v $
 *
 * The class LagrangianTIDS  allows to define  and compute a generic ndof-dimensional
 * Lagrangian Linear Time Invariant Dynamical System of the form :
 * \f[
 * M \ddot q + C \dot q + K q =  F_{Ext}( q , t) + p,
 * \f]
 * where
 *    - \f$q \in R^{ndof} \f$ is the set of the generalized coordinates,
 *    - \f$ \dot q  \in R^{ndof} \f$  the velocity, i. e. the time derivative of the  generalized coordinates.
 *    - \f$ \ddot q  \in R^{ndof} \f$  the acceleration, i. e. the second time derivative of the  generalized coordinates.
 *    - \f$ p  \in R^{ndof} \f$  the forces due to the Non Smooth Interaction. In particular case of Non Smooth evolution,  the variable p stored the impulse and not the force.
 *    -  \f$ M \in  R^{ndof \times ndof} \f$ is Mass matrix stored in the SiconosMatrix mass.
 *    -  \f$ K \in  R^{ndof \times ndof} \f$ is the stiffness matrix stored in the SiconosMatrix K.
 *    -  \f$ C \in  R^{ndof \times ndof} \f$ is the viscosity matrix stored in the SiconosMatrix C.
 *
 *
 *
 *
  *    -  \f$ F_{Int}(\dot q , q , t)\f$ the internal forces stored in the SiconosVector fExt.
 *    -  \f$ F_{Ext}( q , t) \f$ the external forces stored in the SiconosVector fInt.
 *
 *
 * One word on the initial condition.
 *
 * One word on the bilateral constraint
 *
 *
 * The master Class LagrangianNLDS is specified as follows :
 *    -  \f$ M(q) = M q \f$
 *    -  \f$ Q(\dot q, q) = 0 \f$
 *    -  \f$ F_{Int}(\dot q , q , t) = -C \dot q - K q \f$
 *
 *
 *
 * As for the master Class LagrangianNLDS, the state of the master class DynamicalSystem is defined by \f$ x = \left[\begin{array}{c}q \\ \dot q\end{array}\right]\f$ and then \f$ n= 2 ndof \f$ and the VectorField
 * is specified as :
 * \f[
 * f(x,t) = \left[\begin{array}{cc}
 *  0_{ndof \times ndof} & I_{ndof \times ndof} \\
 * M^{-1}(q)\left[   F_{Int}(\dot q , q , t)+F_{Ext}( q , t) -  Q(\dot q, q) \right]\\
 * \end{array}\right]
 * \f]
 *  and the input due to the non smooth law by
 * \f[
 * r = \left[\begin{array}{c}0 \\ p \end{array}\right]
 * \f]
 *
 * \todo Automatically, specify the function of LagrangianNLDS such as
 *          Mass, QNL Inertia , FInt = K q +c velocity,
 */
class LagrangianTIDS : public LagrangianNLDS
{
public:

  /** \fn LagrangianTIDS()
   *  \brief default constructor
   */
  LagrangianTIDS();

  /** \fn LagrangianTIDS(DSXML*)
   *  \brief constructor with XML object of the LagrangianTIDS
   *  \param DSXML* : the XML object corresponding
   */
  LagrangianTIDS(DSXML*);

  ~LagrangianTIDS();

  // getter/setter
  /** \fn SiconosMatrix getK(void)
   *  \brief allow to get the SiconosMatrix K
   *  \return the SiconosMatrix K
   */
  inline SiconosMatrix getK(void) const
  {
    return this->K;
  };

  /** \fn SiconosMatrix getC(void)
   *  \brief allow to get the SiconosMatrix C
   *  \return the SiconosMatrix C
   */
  inline SiconosMatrix getC(void) const
  {
    return this->C;
  };

  /** \fn SiconosMatrix* getKPtr(void)
   *  \brief allow to get the SiconosMatrix K
   *  \return the SiconosMatrix* K
   */
  SiconosMatrix* getKPtr(void);

  /** \fn SiconosMatrix* getCPtr(void)
   *  \brief allow to get the SiconosMatrix C
   *  \return the SiconosMatrix* C
   */
  SiconosMatrix* getCPtr(void);

  /** \fn void setK(SiconosMatrix)
   *  \brief allow to set the SiconosMatrix K
   *  \param the SiconosMatrix to set K
   */
  inline void setK(const SiconosMatrix &K)
  {
    this->K = K;
  };

  /** \fn void setC(SiconosMatrix)
   *  \brief allow to set the SiconosMatrix C
   *  \param the SiconosMatrix to set C
   */
  inline void setC(const SiconosMatrix &C)
  {
    this->C = C;
  };


  //////////////////////////

  /** \fn void init()
   *  \brief initialise value of a Lagrangian TIDS
   */
  virtual void init();

  /** \fn void saveDSToXML()
   *  \brief copy the data of the DS to the XML tree
   *  \exception RuntimeException
   */
  void saveDSToXML();

  /** \fn void display()
   *  \brief print the data to the screen
   */
  void display() const;


  /** \fn void createDynamicalSystem(DSXML * dsXML, int number, int ndof,
        SiconosVector* q0, SiconosVector* velocity0, SiconosMatrix* mass,
        string fExt,SiconosMatrix* K, SiconosMatrix* C)
   *  \brief allows to create the DynamicalSystem with an xml file, or the needed data
   *  \param DSXML * : the XML object for this DynamicalSystem
   *  \param int : the number for this DynamicalSystem
   *  \param int : the dimension of this DynamicalSystem
   *  \param SiconosVector* : the initial coordinates of this DynamicalSystem
   *  \param SiconosVector* : the initial velocity of this DynamicalSystem
   *  \param SiconosMatrix* : the mass of this DynamicalSystem
   *  \param string : the indiaction needed to locate and use the fExt plugin
   *  \param SiconosMatrix* : the matrix K of this DynamicalSystem
   *  \param SiconosMatrix* : the matrix C of this DynamicalSystem
  //  *  \param NSDS * : The NSDS which contains this DynamicalSystem
   *  \exception RuntimeException
   */
  void createDynamicalSystem(DSXML * dsXML, int number = -1, int ndof = -1,
                             SiconosVector* q0 = NULL, SiconosVector* velocity0 = NULL,
                             SiconosMatrix* mass = NULL,
                             string fExt = "BasicPlugin:computeFExt",
                             //        string jacobianQFInt="BasicPlugin:computeJacobianQFInt",
                             //        string jacobianVelocityFInt="BasicPlugin:computeJacobianVelocityFInt",
                             //        string jacobianQQNLInertia="BasicPlugin:computeJacobianQQNLInertia",
                             //        string jacobianVelocityQNLInertia="BasicPlugin:computeJacobianVelocityQNLInertia",
                             //        string QNLlInertia="BasicPlugin:computeQNLInertia",
                             SiconosMatrix* K = NULL, SiconosMatrix* C = NULL); //, NSDS * nsds = NULL);

  static LagrangianTIDS* convert(DynamicalSystem* ds);

protected :
  /** \fn void fillDSWithDSXML()
   *  \brief overload of the function for a LagrangianTIDS
   *  \exception RuntimeException
   */
  void fillDSWithDSXML();


private:
  /** specific matrix for a LagrangianTIDS */
  SiconosMatrix K;

  /** specific matrix for a LagrangianTIDS */
  SiconosMatrix C;

  /*
   * fake of plugin ^^
   */
  friend void LTIDSComputeFExt(int* sizeOfq, double* time, double* qPtr, double* fExt);
};

#endif // LAGRANGIANTIDS_H
//$Log: LagrangianTIDS.h,v $
//Revision 1.32  2005/02/11 17:36:00  charlety
//
//_ little "inspection of code"
//_ basic getters and setters passed inline
//_ getters functions passed const
//
//Revision 1.31  2005/01/31 16:26:20  charlety
//
//_ Added a method named "convert" to classes which inherits from another. This is necessary for Python interface, in order to be able to use down-casting mechanism.
//
//Revision 1.30  2004/09/30 08:35:02  jbarbier
//- fonction of the formalisation : fill..With...XML and link... are now
//"protected" and no more "public"
//
//Revision 1.29  2004/09/22 11:16:28  charlety
//
//_ revision of Doxygen comments in modelformalisation
//
//Revision 1.28  2004/09/22 10:54:43  jbarbier
//- light modification according to the attribute mass of the lagrangian dynamical
//systems. The lagrangianNLDS take always an function from a plugin to compute the
//mass, whereas the lagrangianTIDS needs only a matrix.
//
//- xml input files have been modified in consequence
//
//Revision 1.27  2004/09/10 11:26:14  charlety
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
//Revision 1.26  2004/08/23 14:30:01  jbarbier
//- All the dynamical systems can be created in a comand program and added to a
//NSDS. The save is OK, but the creation of the boundary conditions is not yet
//finished.
//
//Revision 1.25  2004/08/20 15:26:44  jbarbier
//- creation of a Model and save in the XML is ok
//- creation of a NSDS and save in the XML is ok
//- creation of a NonLinearSystemDS and save in the XML is OK
//
//Revision 1.24  2004/08/17 15:12:38  jbarbier
//- methods createDynamicalSystem, createBoundaryCondition, createInteraction,
//createRelation and createNSLaw completed with the required attributes
//
//Revision 1.23  2004/08/12 11:55:14  jbarbier
//- new methods createModel, createNSDS, createStrategy, ...
//they now allow to make the link with upper objects of the platform
//it will be used for the creation of the platform without XML input file
//
//- the createModel method is finished but the attributes of the other objects
//of the platform are missing for the conctruction
//
//Revision 1.22  2004/07/23 14:39:25  jbarbier
//- createModel, createNSDS, createDynamicalSystem, createBoundaryCondition OK
//it's the first step, it do the same thing that before, but the process is
//unified and it must simply add the traitment for the creation of the nodes in
//the DOM tree
//
//Revision 1.21  2004/07/05 12:38:08  charlety
//
//try of false plugin developed in LagrangianTIDS. The Moreau integrator considers it like a LagrangianNLDS, but this is not the plugin which is used to compute the external strength, but a function of the LagrangianTIDS.
//
//Revision 1.20  2004/06/29 08:49:57  acary
//Ajout des commentaires Doxygen et des tages CVS
//