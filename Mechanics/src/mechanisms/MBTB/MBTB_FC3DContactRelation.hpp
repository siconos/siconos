
#ifndef MBTB_FC3DCONTACTRELATION
#define MBTB_FC3DCONTACTRELATION

#include "SiconosKernel.hpp"
#include "MBTB_Contact.hpp"

//! It is a relation dedicated for the unilateral constraint with Coulomb friction.
/*!
  Aggregation to the class MBTB_Contact, the member _pContact contains the CAD informations.
  It derivates from Siconos::NewtonEulerFrom3DLocalFrameR. This class does the link between CAD and Siconos.
 */

class MBTB_FC3DContactRelation : public NewtonEulerFrom3DLocalFrameR
{
protected:
  //! The member containing the CAD properties.
  MBTB_Contact * _pContact;
  MBTB_FC3DContactRelation();
public :
  //! Builder
  /*
    \param [in] a pointer to the MBTB_Contact. Must be allocated/free by the caller.
   */
  MBTB_FC3DContactRelation(MBTB_Contact * _pContact);
  //!This function has to compute the distance between the objects.
  virtual void computeh(double time, Interaction & inter);
  //! Doing nothing.
  virtual ~MBTB_FC3DContactRelation();
};
TYPEDEF_SPTR(MBTB_FC3DContactRelation);

#endif
