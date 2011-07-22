/* Siconos-Kernel, Copyright INRIA 2005-2011.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */

/*! \file SimpleVectorFriends.hpp
  List of friend functions for SimpleVectors.
*/

#ifndef __SimpleVectorFriends__
#define __SimpleVectorFriends__

class VectorNum;

class SimpleMatrix;

/** Copy a subBlock of size sizeB of vIn (from index startIn) into a subBlock
 *  of vOut (from index startOut)
 * \param vIn, a SP::SiconosVector
 * \param vOut, a SP::SiconosVector
 * \param sizeB, an unsigned int
 * \param startIn, an unsigned int
 * \param startOut, an unsigned int
 */
void setBlock(const SiconosVector&, SP::SiconosVector, unsigned int, unsigned int, unsigned int);

/** A==B when (A-B).normInf()<tolerance
 * \param 2 SiconosVector
 * \return a boolean
 */
bool operator ==(const SiconosVector&, const SiconosVector&);

/** multiplication of a vector by a scalar
 *  \param a double
 *  \param a SiconosVector
 *  \return a SimpleVector
 */
SimpleVector operator * (double, const SiconosVector&);

/** multiplication of a vector by a double
 *  \param a SiconosVector
 *  \param a double
 *  \return a SimpleVector
 */
SimpleVector operator * (const SiconosVector&, double);

/** division of the vector by a double
 *  \param a SiconosVector
 *  \param a double
 *  \return a SimpleVector
 *  \exception SiconosVectorException, if the double d = 0
 */
SimpleVector operator / (const SimpleVector&, double);

/** Addition of two vectors
 * \param a SiconosVector
 * \param a SiconosVector
 * \return a SimpleVector
 */
SimpleVector operator + (const SiconosVector&, const SiconosVector&);

/** computes z = x + y
    \param x, a  SiconosVector, IN.
    \param y, a  SiconosVector, IN.
    \param z, a SiconosVector, IN-OUT.
*/
void add(const SiconosVector&, const SiconosVector&, SiconosVector&);

/** Subtraction of two vectors
    \param a SiconosVector (x), IN.
    \param a SiconosVector (y), IN.
    \return a SimpleVector
*/
SimpleVector operator - (const SiconosVector&, const SiconosVector&);

/** computes z = x - y
    \param a SiconosVector (x), IN.
    \param a SiconosVector (y), IN.
    \param a SiconosVector (z), IN-OUT.
*/
void sub(const SiconosVector&, const SiconosVector&, SiconosVector&);

/** computes y = a*x + b*y with atlas axpy.
    \param a, a double.
    \param x, a SiconosVector , IN.
    \param b, a double.
    \param y, a SiconosVector , IN-OUT.
*/
void axpby(double, const SiconosVector&, double, SiconosVector&);

/** computes y = a*x + y with atlas axpy.
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

/** compute the product m1 * trans(m2)
 *  \param 2 SiconosVectors
 *  \return a SimpleMatrix
 */
SimpleMatrix outer_prod(const SiconosVector&, const SiconosVector&);

/** multiplication of a vector by a scalar, y = a*x (init = true) or y += a*x (init = false)
 *  \param a, a double
 *  \param x, a SiconosVector (IN)
 *  \param y, a SiconosVector (IN-OUT)
 *  \param init, a bool, default = true
 */
void scal(double, const SiconosVector&, SiconosVector&, bool = true);

/** multiplication of a vector by a scalar, sub_y = a*sub_x (init = true) or sub_y += a*sub_x (init = false)
 *  \param a, a double
 *  \param x, a SiconosVector (IN)
 *  \param y, a SiconosVector (IN-OUT)
 \param an Index  = [r0x r1x r0y r1y];
 subX is the sub-vector of x, for row numbers between r0x and r1x-1.
 The same for y with riy.
 *  \param init, a bool, default = true
 */
void subscal(double, const SiconosVector&, SiconosVector&, const Index&, bool = true);

/** cross product
 *  \param V1, a SimpleVector of dimention 3.
 *  \param V2, aSimpleVector of dimention 3.
 *  \param VOUT, aSimpleVector of dimention 3, the resulting cross product between V1 and V2.
 */
void cross_product(const SiconosVector&, const SiconosVector&, SiconosVector&);

/** get an absolute vector
 *  \param 1, a SimpleVector (Input).
 *  \param 2, a SimpleVector (Output).
 */

void abs_wise(const SiconosVector&, SiconosVector&);

/** get maximal element of a vector
 *  \param 1, a SimpleVector (Input).
 *  \param 2, a double variable giving the maximum element (Output).
 *  \param 3, an unsigned int variable giving the position of the maximum element (Output)
 */

void getMax(const SiconosVector&, double &, unsigned int &);

/** get minimum element of a vector
 *  \param 1, a SimpleVector (Input).
 *  \param 2, a double variable giving the minimum element (Output).
 *  \param 3, an unsigned int variable giving the position of the minimum element (Output)
 */

void getMin(const SiconosVector&, double &, unsigned int &);


#endif
