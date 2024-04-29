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

#ifndef SubPluggedObject_H
#define SubPluggedObject_H

#include "PluggedObject.hpp"
#include "SiconosException.hpp"
#include "SiconosMatrix.hpp"

namespace siconos::algebra {
// class SiconosMatrix;

}

namespace siconos::plugins {

/*! \file SubPluggedObject.hpp
  \brief utilities for plugin definition to compute extract a column from a matrix computed by
  an user-defined function.
*/

typedef void (*matrixPlugin)(double, unsigned int, unsigned int, double*, unsigned int,
                             double*);

class SubPluggedObject : public PluggedObject {
 private:
  ACCEPT_SERIALIZATION(SubPluggedObject);
  /** function pointer to the parent PluggedObject*/
  void* _parentfPtr{nullptr};
  /** Matrix to temporary storage of the output */
  std::shared_ptr<siconos::algebra::SiconosMatrix> _tmpMat{nullptr};
  /** Column index */
  unsigned int _indx{0};
  /** Number of columns */
  unsigned int _p{0};

  SubPluggedObject() = delete;
  
 public:
  // /** Default Constructor
  //  */
  // SubPluggedObject() = default;

  /** Constructor with the plugin name
   * \param PO a PluggedObject
   * \param n the number of rows of the matrix
   * \param p the number of column of the matrix
   * \param indx the column index (optional)
   */
  SubPluggedObject(const PluggedObject& PO, const unsigned int n, const unsigned int p,
                   const unsigned int indx = 0);

  // /** Copy constructor
  //  * \param SPO a PluggedObject we are going to copy
  //  */
  // SubPluggedObject(const SubPluggedObject& SPO)
  //     : PluggedObject(SPO), _indx(SPO.getIndex()), _p(SPO.getp())
  // {
  //   _parentfPtr = SPO.getParentfPtr();
  //   _tmpMat = std::make_shared<siconos::algebra::SiconosMatrix>(SPO.getTmpMat());
  // }

  /** destructor
   */
  ~SubPluggedObject() noexcept = default;

  /** Intermediate function to compute the column of a plugged matrix
   * \param time current time
   * \param n the length of the vector
   * \param M the vector used for storage
   * \param sizez the size of z
   * \param z a vector used as a parameter
   */
  void computeAndExtract(double time, unsigned int n, double* M, unsigned int sizez,
                         double* z);

  /* Set column index
   * \param newIndx the new column index
   */
  inline void setIndex(unsigned int newIndx) { _indx = newIndx; };

  /* Get column index
   * \return the column index
   */
  inline unsigned int getIndex() const { return _indx; };

  /** Get the number of row
   * \return the number of row
   */
  inline unsigned int getp() const { return _p; };

  /** Get the user defined plugin
   * \return the user defined plugin
   */
  inline void* getParentfPtr() const { return _parentfPtr; };

  /** Get the user defined plugin
   * \return the user defined plugin
   */
  inline const siconos::algebra::SiconosMatrix& getTmpMat() const { return *_tmpMat; };
};
}  // namespace siconos::plugins
#endif
