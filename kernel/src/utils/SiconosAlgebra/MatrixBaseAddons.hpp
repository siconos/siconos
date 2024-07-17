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

/*! \file MatrixBaseAddons.hpp

  Eigen MatrixBase class extensions

  See https://eigen.tuxfamily.org/dox/TopicCustomizing_Plugins.html

*/

// Setters

inline void setValue(Index i, Scalar v) {
  assert(this->cols() == 1);
  this->operator()(i) = v;
}

inline void setValue(Index row, Index col, Scalar v) { this->operator()(row, col) = v; }

/** set a sub-block of the current vector

    \param vectorIndex the beginning of the destination range
    \param vector vector to be copied
*/
template <typename OtherDerived>
inline void setBlock(Index vectorIndex, const MatrixBase<OtherDerived>& vector) {
  assert((vector.cols() == 1) && this->cols() == 1);
  assert(this->rows() >= vectorIndex + vector.rows());
  this->segment(vectorIndex, vector.rows()) = vector;
}

// Set current matrix elements, starting from row row_min and column col_min,
// with the values of the matrix m
template <typename OtherDerived>
inline void setBlock(unsigned int row_min, unsigned int col_min,
                     const MatrixBase<OtherDerived>& m) {
  // assert(m != *this);
  assert(row_min < this->rows() && "row is out of range");
  assert(col_min < this->cols() && "column is out of range");

  unsigned int row_max, col_max;
  row_max = m.rows() + row_min;
  col_max = m.cols() + col_min;
  assert(row_max <= this->rows());
  assert(col_max <= this->cols());

  this->block(row_min, col_min, m.rows(), m.cols()) = m;
}

/** set column col of the current matrix with vector v
 *
 *  \param col index column we want to set
 *  \param vIn a SiconosVector containing the new column
 */
template <typename OtherDerived>
inline void setCol(Index c, const MatrixBase<OtherDerived>& vIn) {
  assert(vIn.cols() == 1 && "setCol(r, vIn) : vIn must have only one column.");
  assert(c < this->cols() && "r is out of range.");
  assert(vIn.size() == this->rows() && "inconsistent sizes between this and v.");

  this->col(c) = vIn;
}

/** set row r of the current matrix with vector v
 *
 *  \param r index row we want to set
 *  \param vIn a SiconosVector containing the new row
 */
template <typename OtherDerived>
inline void setRow(Index r, const MatrixBase<OtherDerived>& vIn) {
  assert(vIn.cols() == 1 && "setCol(r, vIn) : vIn must have only one column.");
  assert(r < this->rows() && "r is out of range.");
  assert(vIn.size() == this->cols() && "inconsistent sizes between this and v.");

  this->row(r) = vIn;
}

// Getters

inline Scalar getValue(Index i, Index j) const {
  assert(this->cols() > 1);
  return this->operator()(i, j);
}
inline Scalar& getValue(Index i, Index j) {
  assert(this->cols() > 1);
  return this->operator()(i, j);
}
inline Scalar getValue(Index i) const {
  assert(this->cols() == 1);
  return this->operator[](i);
}
inline Scalar& getValue(Index i) {
  assert(this->cols() == 1);
  return this->operator[](i);
}

// Misc

inline void eye() { this->setIdentity(); }

inline Scalar normInf() { return this->cwiseAbs().rowwise().sum().maxCoeff(); }
