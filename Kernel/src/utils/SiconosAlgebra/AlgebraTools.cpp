/* Siconos-Kernel, Copyright INRIA 2005-2013.
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

#include "AlgebraTools.hpp"
#include "SiconosMatrix.hpp"
#include <boost/numeric/ublas/io.hpp>
#include "expm.hpp"

namespace Siconos {
  namespace algebra {
    namespace tools {

      void expm(SiconosMatrix& A, SiconosMatrix& Exp, bool computeAndAdd)
      {
        // Implemented only for dense matrices.
        // Note FP : Maybe it works for others but it has not been
        // tested here --> to be done
        // Do not work with sparse.
        assert(Exp.getNum() == 1 || A.getNum() == 1); 
        if(computeAndAdd)
          *Exp.dense() += expm_pad(*A.dense());
        else
          *Exp.dense() = expm_pad(*A.dense());
          
      }
    } // namespace tools
  } // namespace algebra
} // namespace Siconos
