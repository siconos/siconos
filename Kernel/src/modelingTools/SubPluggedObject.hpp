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

#ifndef SubPluggedObject_H
#define SubPluggedObject_H

#include "PluggedObject.hpp"

/*! \file SubPluggedObject.hpp
  \brief utilities for plugin definition to compute extract a column from a matrix computed by an user-defined function.
*/

typedef void (*matrixPlugin)(double, unsigned int, unsigned int, double*, unsigned int, double*);

class SubPluggedObject : public PluggedObject
{
private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(SubPluggedObject);
  /** function pointer to the parent PluggedObject*/
  PluginHandle _parentfPtr;
  /** Matrix to temporary storage of the output */
  SP::SiconosMatrix _tmpMat;
  /** Column index */
  unsigned int _indx;
  /** Number of columns */
  unsigned int _p;

public:

  /** Default Constructor
   */
  SubPluggedObject(): _indx(0), _p(0)
  {
    PluggedObject();
    _parentfPtr = 0;
  };

  /** Constructor with the plugin name
   * \param PO a PluggedObject
   * \param n the number of rows of the matrix
   * \param p the number of column of the matrix
   * \param indx the column index (optional)
   */
  SubPluggedObject(const PluggedObject& PO, const unsigned int n, const unsigned int p, const unsigned int indx = 0):  _indx(indx), _p(p)
  {
    _pluginName = "Sub" + PO.getPluginName();
    _tmpMat.reset(new SimpleMatrix(n, p));
#if (__GNUG__ && !( __clang__ || __INTEL_COMPILER))
    fPtr = (void *)&SubPluggedObject::computeAndExtract;
    _parentfPtr = PO.fPtr;
#else
    RuntimeException::selfThrow("SubPluggedObject must be compiled with GCC !");
#endif
  };

  /** Copy constructor
   * \param SPO a PluggedObject we are going to copy
  */
  SubPluggedObject(const SubPluggedObject& SPO): _indx(SPO.getIndex()), _p(SPO.getp())
  {
    _pluginName = SPO.getPluginName();
    fPtr = SPO.fPtr;
    _parentfPtr = SPO.getParentfPtr();
    _tmpMat.reset(new SimpleMatrix(SPO.getTmpMat()));
  }

  /** destructor
   */
  ~SubPluggedObject() {};

  /** Intermediate function to compute the column of a plugged matrix
   * \param time current time
   * \param n the length of the vector
   * \param M the vector used for storage
   * \param sizez the size of z
   * \param z a vector used as a parameter
   */
  void computeAndExtract(double time, unsigned int n, double* M, unsigned int sizez, double* z)
  {
    ((matrixPlugin)_parentfPtr)(time, n, _p, &(*_tmpMat)(0, 0), sizez, z);
    for (unsigned int i = 0; i < n; i++)
    {
      M[i] = (*_tmpMat)(i, _indx);
    }
  };

  /* Set column index
   * \param newIndx the new column index
   */
  inline void setIndex(unsigned int newIndx)
  {
    _indx = newIndx;
  };

  /* Get column index
   * \return the column index
   */
  inline unsigned int getIndex() const
  {
    return _indx;
  };

  /** Get the number of row
   * \return the number of row
   */
  inline unsigned int getp() const
  {
    return _p;
  };

  /** Get the user defined plugin
   * \return the user defined plugin
   */
  inline PluginHandle getParentfPtr() const
  {
    return _parentfPtr;
  };

  /** Get the user defined plugin
   * \return the user defined plugin
   */
  inline const SiconosMatrix& getTmpMat() const
  {
    return *_tmpMat;
  };

};
TYPEDEF_SPTR(SubPluggedObject);
#endif
