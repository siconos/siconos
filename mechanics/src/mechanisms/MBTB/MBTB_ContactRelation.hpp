#ifndef MBTB_CONTACTRELATION
#define MBTB_CONTACTRELATION

#include "MechanicsFwd.hpp"
#include <SiconosKernel.hpp>
#include "MBTB_Contact.hpp"
//! It is a relation dedicated for the simple unilateral (ie: without  Coulomb friction).
/*!
  Aggregation to the class MBTB_Contact, the member _pContact contains the CAD informations.
  It derivates from Siconos::NewtonEulerFrom1DLocalFrameR. This class does the link between CAD and Siconos.
 */
class MBTB_ContactRelation : public NewtonEulerFrom1DLocalFrameR
{

protected:
  MBTB_Contact * _pContact;
  MBTB_ContactRelation();
public:
  /** Constructor
   *  \param pC [in] a pointer to the MBTB_Contact. Must be allocated/free by the caller.
   */
  MBTB_ContactRelation(MBTB_Contact * pC);
  /** This function has to compute the distance between the objects.
   * \param time  the given time
   * \param q0 the position
   * \param y the output
   */
  virtual void computeh(double time, BlockVector& q0, SiconosVector& y);
  //! Doing nothing.
  virtual ~MBTB_ContactRelation();

  ACCEPT_STD_VISITORS();

};


#endif
