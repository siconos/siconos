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


#include "SiconosException.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosMatrixOp.hpp"  // For setBlock


void siconos::algebra::setBlock(const SiconosMatrix &input_matrix,
                                std::shared_ptr<SiconosMatrix> output_matrix,
                                const std::vector<std::size_t> &dim,
                                const std::vector<std::size_t> &start) {
  // To copy a subBlock of input_matrix into a subBlock of output_matrix.
  // dim[0], dim[1]: number of rows and columns of the sub-block
  // start[0], start[1]: position (row, column) of the first element of the subBlock in
  // input_matrix start[2], start[3]: position (row, column) of the first element of the
  // subBlock in output_matrix

  // if (input_matrix == output_matrix) // useless op => nothing to be done.
  //   return;
    
    // if(output_matrix->isZero() || output_matrix->isIdentity())
    //     THROW_EXCEPTION("output_matrix is read-only (zero or identity matrix?)."); TODO : deal with this, as there are no read-only matrices in eigen.

    // Check dimension
    std::vector<std::size_t> MDim(4);  // dim. of matrices input_matrix and output_matrix.
    MDim[0] = input_matrix.size(0);
    MDim[1] = input_matrix.size(1);
    MDim[2] = output_matrix->size(0);
    MDim[3] = output_matrix->size(1);

    for (unsigned int i = 0; i < 4; ++i)
        if (start[i] >= MDim[i])
        THROW_EXCEPTION(
            "matrices, setBlock(input_matrix, ...): sub-block indices are out of range.");

    // index position of the last element in subBlock ...
    std::vector<std::size_t> end(4);
    end[0] = dim[0] + start[0];
    end[1] = dim[1] + start[1];
    end[2] = dim[0] + start[2];
    end[3] = dim[1] + start[3];

    for (unsigned int i = 0; i < 4; ++i)
        if (end[i] > MDim[i]) THROW_EXCEPTION("sub-block last indices are out of range.");

    // Elements from row/col start[i] to row/col (end[i]-1) will be copied.

    output_matrix->block(start[2], start[3], dim[0], dim[1]) = input_matrix.block(start[0], start[1], dim[0], dim[1]);
}


// void siconos::algebra::setBlock(const SiconosMatrix &input_matrix,
//                                 std::shared_ptr<MapType> output_matrix,
//                                 const std::vector<std::size_t> &dim,
//                                 const std::vector<std::size_t> &start) {
//   // To copy a subBlock of input_matrix into a subBlock of output_matrix.
//   // dim[0], dim[1]: number of rows and columns of the sub-block
//   // start[0], start[1]: position (row, column) of the first element of the subBlock in
//   // input_matrix start[2], start[3]: position (row, column) of the first element of the
//   // subBlock in output_matrix

//   // if (input_matrix == output_matrix) // useless op => nothing to be done.
//   //   return;
    
//     // if(output_matrix->isZero() || output_matrix->isIdentity())
//     //     THROW_EXCEPTION("output_matrix is read-only (zero or identity matrix?)."); TODO : deal with this, as there are no read-only matrices in eigen.

//     // Check dimension
//     std::vector<std::size_t> MDim(4);  // dim. of matrices input_matrix and output_matrix.
//     MDim[0] = input_matrix.size(0);
//     MDim[1] = input_matrix.size(1);
//     MDim[2] = output_matrix->rows();
//     MDim[3] = output_matrix->cols();

//     for (unsigned int i = 0; i < 4; ++i)
//         if (start[i] >= MDim[i])
//         THROW_EXCEPTION(
//             "matrices, setBlock(input_matrix, ...): sub-block indices are out of range.");

//     // index position of the last element in subBlock ...
//     std::vector<std::size_t> end(4);
//     end[0] = dim[0] + start[0];
//     end[1] = dim[1] + start[1];
//     end[2] = dim[0] + start[2];
//     end[3] = dim[1] + start[3];

//     for (unsigned int i = 0; i < 4; ++i)
//         if (end[i] > MDim[i]) THROW_EXCEPTION("sub-block last indices are out of range.");

//     // Elements from row/col start[i] to row/col (end[i]-1) will be copied.

//     output_matrix->block(start[2], start[3], dim[0], dim[1]) = input_matrix.block(start[0], start[1], dim[0], dim[1]);
// }
