
#ifndef MBTB_FC3DCONTACTRELATION
#define MBTB_FC3DCONTACTRELATION

#include "MechanismsFwd.hpp"
#include <SiconosKernel.hpp>
#include "MBTB_Contact.hpp"

//! It is a relation dedicated for the unilateral constraint with Coulomb friction.
/*!
  Aggregation to the class MBTB_Contact, the member _pContact contains the CAD information.
  It derivates from Siconos::NewtonEuler3DR. This class does the link between CAD and Siconos.
 */

class MBTB_FC3DContactRelation : public NewtonEuler3DR
{
protected:
  //! The member containing the CAD properties.
  MBTB_Contact * _pContact;
  MBTB_FC3DContactRelation();
public :
  /** Constructor
   * \param _pContact [in] a pointer to the MBTB_Contact. Must be allocated/free by the caller.
   */
  MBTB_FC3DContactRelation(MBTB_Contact * _pContact);
  
  /** This function has to compute the distance between the objects.
   * \param time the given  time
   * \param q0 the position
   * \param y the output
   */
  virtual void computeh(double time, const BlockVector& q0, SiconosVector& y);
  
  //! Doing nothing.
  virtual ~MBTB_FC3DContactRelation();

  ACCEPT_STD_VISITORS();

};

#endif
