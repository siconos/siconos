/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef MBTB_CONTACT
#define MBTB_CONTACT

#include <memory>
#include <string>

namespace siconos::modeling {

class NewtonEuler1DR;
class Interaction;
}  // namespace siconos::modeling

namespace siconos::mechanisms {

/**
   Contact class.
   This class describes a CAD contact: it contains the identifier of the CAD
   model of the contact. It is either only unilateral or unilateral with
   Coulomb friction, depending on withFriction parameter of the builder. It
   builds the corresponding member _Relation either MBTB_ContactRelation or
   MBTB_FC3DContactRelation.
*/
class MBTB_Contact {
  //! For the unilateral case.
  friend class MBTB_ContactRelation;
  //! For the Coulomb friction case.
  friend class MBTB_FC3DContactRelation;

 protected:
  //! The contact name.
  std::string _ContactName{""};

  double _dist{0.};

  /*!
    To avoid unuseful call at the same date.
   */
  double _curTimeh{-1.};
  /*!
    Built during the building of the MBTB_Contact.
    A link to the relation, either MBTB_ContactRelation or
    MBTB_FC3DContactRelation.
   */
  std::shared_ptr<siconos::modeling::NewtonEuler1DR> _Relation{nullptr};
  /*!
    Built during the building of the MBTB_Contact.
    A link to the interaction.
   */
  std::shared_ptr<siconos::modeling::Interaction> _interaction{nullptr};

  //! The id of the contact. -1 means unset
  int _id{-1};

  //! The index of the body carrying the first face of the contact.
  int _indexBody1{-1};

  //! The index of the body carrying the second face of the contact. -1 means
  //! the face doesn't move.
  int _indexBody2{-1};

  //! The index of the cad model.
  int _indexCAD1{-1};

  //! The index of the cad model.
  int _indexCAD2{-1};

  //! The value 0 means without friction 1 means with.
  bool _withFriction{false};

  //! The normal restitution coefficient.
  double _en{0.5};

  //! The tangential restitution coefficient(not used).
  double _et{0.};

  //! The offset sub to the distance computation.
  double _Offset{0.01};

  //! If normalFromFace1, the normal is computed from the face1 else from the
  //! face2.
  bool _normalFromFace1{true};

  //! To know if P1 is translated of _Offset*N or P2.
  bool _OffsetP1{true};

  // Rule of five
  MBTB_Contact() = delete;
  MBTB_Contact(const MBTB_Contact &) = delete;
  MBTB_Contact(MBTB_Contact &&) = delete;
  MBTB_Contact &operator=(const MBTB_Contact &) = delete;
  MBTB_Contact &operator=(MBTB_Contact &&) = delete;

 public:
  /** Builds the member _Relation either MBTB_ContactRelation or

      \param[in] id the identier of the contact
      \param[in] contactName name of the contact
      \param[in] indexBody1 index of the first body
      \param[in] indexBody2 index of the second body
         If negatif, it means the body is not movable (ex: the ground)
      \param[in] indexCAD1 the index of the cad model
         used to define the surface of the indexBody1
      \param[in] indexCAD2 the index of the cad model
         used to define the surface of the indexBody2
      \param[in] withFriction if 0, the contact is without friction.
   */
  MBTB_Contact(int id, const std::string &contactName, int indexBody1,
               int indexBody2, int indexCAD1, int indexCAD2, bool withFriction);

  virtual ~MBTB_Contact() noexcept = default;

  /** \return normal restitution coefficient */
  auto en() const { return _en; };

  /** \return tangential restitution coefficient */
  auto et() const { return _et; };

  /** set normal restitution coefficient */
  void set_en(double en) { _en = en; };

  /** set tangential restitution coefficient */
  void set_et(double et) { _et = et; };

  /** \return index of the body carrying the first face of the contact*/
  auto indexBody1() const { return _indexBody1; }

  /** \return index of the body carrying the second face of the contact*/
  auto indexBody2() const { return _indexBody2; }

  /** \return index of the CAD model */
  auto indexCAD1() const { return _indexCAD1; }

  /** \return index of the CAD model */
  auto indexCAD2() const { return _indexCAD2; }

  /** Set if offset must be applied to P1  (or not) */
  void set_offset_to_P1(bool val) { _OffsetP1 = val; }

  /** \return True if an offset must be applied to P1 */
  auto offset_to_P1() const { return _OffsetP1; }

  /** Set offset value */
  void set_offset(double val) { _Offset = val; }

  /** \return offset value */
  auto offset() const {return _Offset;}


  /** Set to true to compute the normal from face1 */
  void set_normal_from_face1(bool val){_normalFromFace1 = val;}
  
  /** \return a pointer to the relation
   */
  inline auto relation() { return _Relation; }

  inline std::shared_ptr<siconos::modeling::Interaction> interaction() {
    return _interaction;
  }
  void setInteraction(
      std::shared_ptr<siconos::modeling::Interaction> newInteraction);

  /** To get the name of the contact.
   * \return char * contactName
   */
  inline auto contactName() { return _ContactName; }
};
}  // namespace siconos::mechanisms
#endif
