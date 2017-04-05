/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
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

/*! \file SiconosVectorFriends.hpp
  List of friend functions for SiconosVectors.
*/

#ifndef __SiconosVectorFriends__
#define __SiconosVectorFriends__


/** Copy a subBlock of size sizeB of vIn (from index startIn) into a subBlock
 *  of vOut (from index startOut)
 * \param vIn block to copy
 * \param vOut vector to change (destination)
 * \param sizeB size of the block to copy
 * \param startIn starting position for the block (vIn)
 * \param startOut starting position for the destination (vOut)
 */
void setBlock(const SiconosVector& vIn, SP::SiconosVector vOut, unsigned int sizeB, unsigned int startIn, unsigned int startOut);

/** A==B when (A-B).normInf()<tolerance
 * \param 2 SiconosVector
 * \return a boolean
 */
bool operator ==(const SiconosVector&, const SiconosVector&);

/** multiplication of a vector by a scalar
 *  \param a double
 *  \param a SiconosVector
 *  \return a SiconosVector
 */
SiconosVector operator * (double, const SiconosVector&);

/** multiplication of a vector by a double
 *  \param a SiconosVector
 *  \param a double
 *  \return a SiconosVector
 */
SiconosVector operator * (const SiconosVector&, double);

/** division of the vector by a double
 *  \param a SiconosVector
 *  \param a double
 *  \return a SiconosVector
 *  \exception SiconosVectorException, if the double d = 0
 */
SiconosVector operator / (const SiconosVector&, double);

/** Addition of two vectors
 * \param a SiconosVector
 * \param a SiconosVector
 * \return a SiconosVector
 */
SiconosVector operator + (const SiconosVector&, const SiconosVector&);

/** computes z = x + y
    \param x, a  SiconosVector, IN.
    \param y, a  SiconosVector, IN.
    \param z, a SiconosVector, IN-OUT.
*/
void add(const SiconosVector&, const SiconosVector&, SiconosVector&);

/** Subtraction of two vectors
    \param a SiconosVector (x), IN.
    \param a SiconosVector (y), IN.
    \return a SiconosVector
*/
SiconosVector operator - (const SiconosVector&, const SiconosVector&);

/** computes z = x - y
    \param a SiconosVector (x), IN.
    \param a SiconosVector (y), IN.
    \param a SiconosVector (z), IN-OUT.
*/
void sub(const SiconosVector&, const SiconosVector&, SiconosVector&);

/** computes y = a*x + b*y with blas axpy.
    \param a, a double.
    \param x, a SiconosVector , IN.
    \param b, a double.
    \param y, a SiconosVector , IN-OUT.
*/
void axpby(double, const SiconosVector&, double, SiconosVector&);

/** computes y = a*x + y with blas axpy.
    \param a, a double.
    \param x, a SiconosVector , IN.
    \param y, a SiconosVector , IN-OUT.
*/
void axpy(double, const SiconosVector&, SiconosVector&);

/** compute dot product m1.m2
 *  \param 2 SiconosVectors
 *  \return a double
 */
double inner_prod(const SiconosVector&, const SiconosVector&);

///** compute the product m1 * trans(m2)
// *  \param 2 SiconosVectors
// *  \return a SimpleMatrix
// */
//SimpleMatrix outer_prod(const SiconosVector&, const SiconosVector&);

/** multiplication of a vector by a scalar, y = a*x (init = true) or y += a*x (init = false)
 *  \param a a double
 *  \param[in] x a SiconosVector
 *  \param[in,out] y a SiconosVector
 *  \param init if true y = a*x else y += a*x (default = true)
 */
void scal(double a, const SiconosVector& x, SiconosVector& y, bool init = true);

/** multiplication of a vector by a scalar, sub_y = a*sub_x (init = true) or sub_y += a*sub_x (init = false)
 * \param a a double
 * \param x a SiconosVector (IN)
 * \param y a SiconosVector (IN-OUT)
 * \param coord an Index  = [r0x r1x r0y r1y];
 subX is the sub-vector of x, for row numbers between r0x and r1x-1.
 The same for y with riy.
 * \param init if true sub_y = a*sub_x else sub_y += a*sub_x (default true)
 */
void subscal(double a, const SiconosVector& x, SiconosVector& y, const Index& coord, bool init = true);

/** cross product
 *  \param V1 a SiconosVector of dimention 3.
 *  \param V2 aSiconosVector of dimention 3.
 *  \param VOUT aSiconosVector of dimention 3, the resulting cross product between V1 and V2.
 */
void cross_product(const SiconosVector& V1, const SiconosVector& V2, SiconosVector& VOUT);

/** get an absolute vector
 *  \param V a SiconosVector (Input).
 *  \param Vabs a SiconosVector (Output)
 */

void abs_wise(const SiconosVector& V, SiconosVector& Vabs);

/** get maximal element of a vector
 *  \param 1, a SiconosVector (Input).
 *  \param 2, a double variable giving the maximum element (Output).
 *  \param 3, an unsigned int variable giving the position of the maximum element (Output)
 */

void getMax(const SiconosVector&, double &, unsigned int &);

/** get minimum element of a vector
 *  \param 1, a SiconosVector (Input).
 *  \param 2, a double variable giving the minimum element (Output).
 *  \param 3, an unsigned int variable giving the position of the minimum element (Output)
 */

void getMin(const SiconosVector&, double &, unsigned int &);


#endif
