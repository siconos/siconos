/* Siconos-IO, Copyright INRIA 2005-2011.
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

/*! \file Dump.hpp
  \brief provides pre-compiled functions for a full Siconos Model
  serialization through libSiconosIO library */

#ifndef Dump_hpp
#define Dump_hpp

#include "SiconosFull.hpp"


/** \namespace Siconos::IO
    the place for non templated Siconos IO functions
*/
namespace Siconos
{
namespace IO
{

/** save a Siconos Model with the full simulation state into a
 *  file
 * \param model
 * \param filename with extension : .xml, .dat (binary archive)
 */
void save(SP::Model model, std::string filename);

/** load a Siconos Model with the full simulation state from file
 * \param filename
 * \return a SP::Model
 */
SP::Model load(std::string filename);

}
}


#endif
