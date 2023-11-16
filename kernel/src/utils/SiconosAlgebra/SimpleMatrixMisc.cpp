/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
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

// #include <boost/numeric/bindings/lapack.hpp>
// #include <boost/numeric/bindings/std/vector.hpp>
// #include <boost/numeric/bindings/ublas/matrix.hpp>
// #include <boost/numeric/bindings/ublas/symmetric.hpp>
// #include <boost/numeric/bindings/ublas/vector.hpp>
// #include <boost/numeric/ublas/banded.hpp>
// #include <boost/numeric/ublas/io.hpp>
// #include <boost/numeric/ublas/lu.hpp>
// #include <boost/numeric/ublas/matrix.hpp>
// #include <boost/numeric/ublas/matrix_proxy.hpp>
// #include <boost/numeric/ublas/matrix_sparse.hpp>
// #include <boost/numeric/ublas/symmetric.hpp>
// #include <boost/numeric/ublas/triangular.hpp>
// #include <boost/numeric/ublas/vector.hpp>
// #include <boost/numeric/ublas/vector_proxy.hpp>

#include "BlockMatrix.hpp"
#include "SiconosException.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosMatrixOp.hpp"
#include "SiconosMatrixVectorOp.hpp"
#include "SiconosVector.hpp"
#include "SimpleMatrix.hpp"
#include "Simulation.hpp"
#include "Tools.hpp"  // enum_to_string
#include "determinant.hpp"
#include "expm.hpp"  // boost contribs expm_pad
#include "io.hpp"


void siconos::algebra::normInfByColumn(const SimpleMatrix &m, SiconosVector &v)
{
    if(v.size() != m.size(1)) THROW_EXCEPTION("the given vector does not have the right length");
    for (size_t i = 0; i < m.size(1); i++)
    {
        v(i) = m.col(i).norm();
    }
}

bool siconos::algebra::checkSymmetry(SiconosMatrix &m, double tol) {
    auto m_trans = m.transpose();
    auto tmp = m - m_trans;
    double err = tmp.normInf();
    if (m_trans.normInf() > 0.0) {
    err /= m_trans.normInf();
    }
    // std::cout << "err_rel  ="<< err <<"\n";
    return (err < tol);
}


siconos::algebra::SimpleMatrix siconos::algebra::readMatrixFromFile(const std::string &filename, bool ascii) {
    SimpleMatrix m;
    if (ascii) {
        io::read(filename, m, io::ASCII_IN);
    } else {
        io::read(filename, m, io::BINARY_IN);
    }
    return m;
}