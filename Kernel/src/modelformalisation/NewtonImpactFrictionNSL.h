//$Id: NewtonImpactFrictionNSL.h,v 1.1 2005/03/22 15:55:05 jbarbier Exp $
#ifndef NEWTONIMPACTFRICTIONNSLAW_H
#define NEWTONIMPACTFRICTIONNSLAW_H

#include "NonSmoothLaw.h"
#include "NewtonImpactFrictionNSLXML.h"

/** \class NewtonImpactFrictionNSL
 *  \brief Specific NonSmoothLaw for the Newton impact friction model
 *  \author Jean-Michel Barbier
 *  \version 0.1
 *  \date (Creation) March 22, 2005
 *
 * $Date: 2005/03/22 15:55:05 $
 * $Revision: 1.1 $
 * $Author: jbarbier $
 * $Source: /CVS/Siconos/SICONOS/src/modelformalisation/NewtonImpactFrictionNSL.h,v $
 *
 */


class NewtonImpactFrictionNSL : public NonSmoothLaw
{

public:

  /** \fn NewtonImpactFrictionNSL()
   *  \brief default constructor
   */
  NewtonImpactFrictionNSL();

  /** \fn NewtonImpactFrictionNSL(NonSmoothLawXML*)
   *  \brief constructor with XML object of the NewtonImpactFrictionNSL
   *  \param NonSmoothLawXML* : the XML object corresponding
   */
  NewtonImpactFrictionNSL(NonSmoothLawXML*);

  /** \fn NewtonImpactFrictionNSL(double en, double et, double mu)
   *  \brief constructor with the value of the NewtonImpactFrictionNSL attributes
   *  \param double : normal e coefficient
   *  \param double : tangent e coefficient
   *  \param double : friction coefficient
   */
  NewtonImpactFrictionNSL(double en, double et, double mu);

  ~NewtonImpactFrictionNSL();

  /** \fn bool isVerified(void);
   *  \brief check the ns law to see if it is verified
   *  \return a boolean value whioch determines if the NS Law is verified
   */
  bool isVerified(void) const;

  /** \fn double getEn(void)
   *  \brief getter of en
   *  \return the value of en
   */
  inline double getEn(void) const
  {
    return this->en;
  };

  /** \fn void setEn(double)
   *  \brief setter of en
   *  \param a double to set en
   */
  inline void setEn(const double En)
  {
    this->en = En;
  };

  /** \fn double getEt(void)
   *  \brief getter of et
   *  \return the value of et
   */
  inline double getEt(void) const
  {
    return this->et;
  };

  /** \fn void setEt(double)
   *  \brief setter of et
   *  \param a double to set et
   */
  inline void setEt(const double Et)
  {
    this->en = Et;
  };

  /** \fn double getMu(void)
   *  \brief getter of mu
   *  \return the value of mu
   */
  inline double getMu(void) const
  {
    return this->mu;
  };

  /** \fn void setMu(double)
   *  \brief setter of mu
   *  \param a double to set mu
   */
  inline void setMu(const double Mu)
  {
    this->mu = Mu;
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

  /** \fn void createNonSmoothLaw(NewtonImpactFrictionNSLXML * nslawXML, double en, double et, double mu)
   *  \brief allows to create the NSLaw with an xml file, or the needed data
   *  \param NewtonImpactLawNSLXML * : the XML object for this NSLaw
   *  \param double : the value of e for this NSLaw
   *  \exception RuntimeException
   */
  void createNonSmoothLaw(NewtonImpactFrictionNSLXML * nslawXML, double en = -1, double et = -1.0, double mu = -1.0);


  /** \fn NewtonImpactFrictionNSL* convert (NonSmoothLaw* nsl)
   *  \brief encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param NonSmoothLaw* : the law which must be converted
   * \return a pointer on the law if it is of the right type, NULL otherwise
   */
  static NewtonImpactFrictionNSL* convert(NonSmoothLaw* nsl);

protected:
  /** \fn void fillNonSmoothLawWithNonSmoothLawXML()
   *  \brief uses the NewtonImpactLawNSLXML of the ComplementarityConditionNSL to fill the fields of this ComplementarityConditionNSL
   *  \exception RuntimeException
   */
  void fillNonSmoothLawWithNonSmoothLawXML();


private:
  /**  \brief The Newton coefficient of restitution
   */
  double en;
  double et;
  double mu;
};

#endif // NewtonImpactFrictionNSL_H
//$Log: NewtonImpactFrictionNSL.h,v $
//Revision 1.1  2005/03/22 15:55:05  jbarbier
//- class NewtonImpactFriction non smooth law added to the kernel
//
//- xml schema modified for this new class
//- xml schema modified to accept a "joker" for further use of a LMGC90 mechanical plugin
//
//- new test added for the loading/saving of a NewtonImpactFrictionNSL
//
