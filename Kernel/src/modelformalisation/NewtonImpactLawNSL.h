//$Id: NewtonImpactLawNSL.h,v 1.10 2005/02/11 17:36:01 charlety Exp $
#ifndef NEWTONIMPACTLAWNSLAW_H
#define NEWTONIMPACTLAWNSLAW_H

#include "NonSmoothLaw.h"
#include "NewtonImpactLawNSLXML.h"

/** \class NewtonImpactLawNSL
 *  \brief Specific NonSmoothLaw for the Newton impact model
 *  \author V. ACARY
 *  \version 0.1
 *  \date (Creation) June 29, 2004
 *
 * $Date: 2005/02/11 17:36:01 $
 * $Revision: 1.10 $
 * $Author: charlety $
 * $Source: /CVS/Siconos/SICONOS/src/modelformalisation/NewtonImpactLawNSL.h,v $
 *
 * The class formalizes the Newton Impact law together with a complementarity condition. i.e.
 * \f[
 * \left\{\begin{array}{l}
 * y \geq 0, \lambda \geq 0, y^{T} \lambda=0\\
 *  if y \leq 0 \quad \mbox{then} \quad \dot y(t^{+}) - e \dot y(t^{-}) \geq 0,   \lambda \geq 0, (\dot y(t^{+}) - e \dot y(t^{-}))^{T} \lambda=0
 * \end{array}\right.
 * \f]
 *
 *
 */


class NewtonImpactLawNSL : public NonSmoothLaw
{

public:

  /** \fn NewtonImpactLawNSL()
   *  \brief default constructor
   */
  NewtonImpactLawNSL();

  /** \fn NewtonImpactLawNSL(NonSmoothLawXML*)
   *  \brief constructor with XML object of the NewtonImpactLawNSL
   *  \param NonSmoothLawXML* : the XML object corresponding
   */
  NewtonImpactLawNSL(NonSmoothLawXML*);

  /** \fn NewtonImpactLawNSL(double c, double d)
   *  \brief constructor with the value of the NewtonImpactLawNSL attributes
   *  \param a double value c
   *  \param a double value d
   */
  NewtonImpactLawNSL(double e);

  ~NewtonImpactLawNSL();

  /** \fn bool isVerified(void);
   *  \brief check the ns law to see if it is verified
   *  \return a boolean value whioch determines if the NS Law is verified
   */
  bool isVerified(void) const;

  /** \fn double getE(void)
   *  \brief getter of e
   *  \return the value of e
   */
  inline double getE(void) const
  {
    return this->e;
  };

  /** \fn void setE(double)
   *  \brief setter of e
   *  \param a double to set e
   */
  inline void setE(const double E)
  {
    this->e = E;
  };

  //////////////////////

  /** \fn void saveNonSmoothLawToXML()
   *  \brief copy the data of the NonSmoothLaw in the XML tree
   */
  void saveNonSmoothLawToXML();

  /** \fn void display()
   *  \brief print the data to the screen
   */
  void display() const;

  /** \fn void createNonSmoothLaw(NewtonImpactLawNSLXML * nslawXML, double e)
   *  \brief allows to create the NSLaw with an xml file, or the needed data
   *  \param NewtonImpactLawNSLXML * : the XML object for this NSLaw
   *  \param double : the value of e for this NSLaw
   *  \exception RuntimeException
   */
  void createNonSmoothLaw(NewtonImpactLawNSLXML * nslawXML, double e = -1); //, Interaction * interaction = NULL);


  /** \fn NewtonImpactLawNSL* convert (NonSmoothLaw* nsl)
   *  \brief encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param NonSmoothLaw* : the law which must be converted
   * \return a pointer on the law if it is of the right type, NULL otherwise
   */
  static NewtonImpactLawNSL* convert(NonSmoothLaw* nsl);

protected:
  /** \fn void fillNonSmoothLawWithNonSmoothLawXML()
   *  \brief uses the NewtonImpactLawNSLXML of the ComplementarityConditionNSL to fill the fields of this ComplementarityConditionNSL
   *  \exception RuntimeException
   */
  void fillNonSmoothLawWithNonSmoothLawXML();


private:
  /**  \brief The Newton coefficient of restitution
   */
  double e;
};

#endif // NewtonImpactLawNSL_H
//$Log: NewtonImpactLawNSL.h,v $
//Revision 1.10  2005/02/11 17:36:01  charlety
//
//_ little "inspection of code"
//_ basic getters and setters passed inline
//_ getters functions passed const
//
//Revision 1.9  2005/01/31 16:26:22  charlety
//
//_ Added a method named "convert" to classes which inherits from another. This is necessary for Python interface, in order to be able to use down-casting mechanism.
//
//Revision 1.8  2004/09/30 08:35:02  jbarbier
//- fonction of the formalisation : fill..With...XML and link... are now
//"protected" and no more "public"
//
//Revision 1.7  2004/09/22 11:16:28  charlety
//
//_ revision of Doxygen comments in modelformalisation
//
//Revision 1.6  2004/09/03 14:41:47  jbarbier
//- new functions to create the boundary condition of the dynamical systems
//- new functions to add an interaction to a NSDS
//- new functions to create the relation and the non-smooth law of an interaciton
//
//Revision 1.5  2004/08/17 15:12:43  jbarbier
//- methods createDynamicalSystem, createBoundaryCondition, createInteraction,
//createRelation and createNSLaw completed with the required attributes
//
//Revision 1.4  2004/08/12 11:55:17  jbarbier
//- new methods createModel, createNSDS, createStrategy, ...
//they now allow to make the link with upper objects of the platform
//it will be used for the creation of the platform without XML input file
//
//- the createModel method is finished but the attributes of the other objects
//of the platform are missing for the conctruction
//
//Revision 1.3  2004/08/11 14:43:45  jbarbier
//- beginning of the mechanism of creation without XML input file of the objects of the platform with the
//creatObjtect methods
//
//- function saveWToXML for Moreau integrator, and same specific functions to save
//M,q and Q,p for LCP and QP
//
//- function to check coherency of the Model
//
//Revision 1.2  2004/07/27 09:32:43  jbarbier
//- functions createNSLaw for complementarityConditionNSL, newtonImpactLawNSL and RelayNSL
//
//Revision 1.1  2004/07/06 14:54:48  acary
//Renaming NSLaw into NonSmoothLaw
//Renaming RelayNSLaw into RelayNSL
//Renaming CCNSLaw into ComplementarityConditionNSL
//Renaming NewtonImpactLaw into NewtonImpactLawNSL
//
//Revision 1.2  2004/07/02 14:48:29  acary
//Added MACRO IN and OUT
//
//Revision 1.1  2004/06/30 09:44:35  acary
//Added NewtonImpactLawNSL
//