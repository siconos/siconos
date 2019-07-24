#ifndef MBTB_CONTACT
#define MBTB_CONTACT
#include "SiconosKernel.hpp"
class MBTB_ContactRelation;
class MBTB_FC3DContactRelation;
//! A Contact class.
/*!
  This class describes a CAD contact: it contains the identifier of the CAD model of the contact.
  It is either only unilateral or unilateral with Coulomb friction, depending on withFriction parameter of the builder.
  It builds the corresponding member _Relation either MBTB_ContactRelation or MBTB_FC3DContactRelation.
 */
class MBTB_Contact
{
  //! For the unilateral case.
  friend class MBTB_ContactRelation;
  //! For the Coulomb friction case.
  friend class MBTB_FC3DContactRelation;
protected:
  //! The contact name.
  char  _ContactName[256];

  double _dist;
  /*!
    To avoid unuseful call at the same date.
   */
  double _curTimeh;
  /*!
    Built during the building of the MBTB_Contact.
    A link to the relation, either MBTB_ContactRelation or MBTB_FC3DContactRelation.
   */
  SP::NewtonEuler1DR _Relation;
  /*!
    Built during the building of the MBTB_Contact.
    A link to the interaction.
   */
  SP::Interaction _interaction;
 
  MBTB_Contact();
public:

  //!Builder. It builds the member _Relation either MBTB_ContactRelation or MBTB_FC3DContactRelation.
  /*!
    \param [in] id an unsigned int, the identier of the contact.
    \param [in] contactName the name of the contact, is copied in the member _ContactName.
    \param [in] indexBody1 an unsigned int, index of the first body carrying the surface of contact.
    \param [in] indexBody2 an int, index of the first body carrying the surface of contact. If negatif, it means the body is the not movable (ex: the ground).
    \param [in] indexCAD1 an unsigned int the index of the cad model defined the surface of the indexBody1.
    \param [in] indexCAD2 an unsigned int the index of the cad model defined the surface of the indexBody2.
    \param [in] withFriction an int, if 0, the contact is without friction.
   */
  MBTB_Contact(unsigned int id,const std::string& contactName, unsigned int indexBody1, int indexBody2,unsigned int indexCAD1,unsigned int indexCAD2,int withFriction);

  /** To get the relation. 
   * \return SP::NewtonEuler1DR the SP on the relation
   */
  inline SP::NewtonEuler1DR relation()
  {
    return _Relation;
  }
  inline SP::Interaction interaction()
  {
    return _interaction;
  }
  void setInteraction(SP::Interaction newInteraction);

  //! The id of the contact.
  unsigned int _id;

  //! The index of the body carrying the first face of the contact.
  unsigned int _indexBody1;

  //! The index of the body carrying the second face of the contact. -1 means the face doesn't move.
  int _indexBody2;

  //! The index of the cad model.
  unsigned int _indexCAD1;

  //! The index of the cad model.
  unsigned int _indexCAD2;

  //! The value 0 means without friction 1 means with.
  int _withFriction;

  //! The normal restitution coefficient.
  double _en;

  //!The tangential restitution coefficient(not used).
  double _et;

  //!The offset sub to the distance computation.
  double _Offset;

  //!If normalFromFace1, the normal is computed from the face1 else from the face2.
  unsigned int _normalFromFace1;

  //! To know if P1 is trnaslated of _Offset*N or P2.
  unsigned int _OffsetP1;

  /** To get the name of the contact.
   * \return char * contactName
   */
  inline char * contactName()
  {
    return _ContactName;
  }
};
#endif
