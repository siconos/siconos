/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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

/*
  Tests functions of Jordan algebra.

 */

#include "JordanAlgebra.h"

static int Arrow_repr_3d_test()
{
    int varsCount = 2;
    int dimension = 3;
    int vecSize = varsCount * dimension;

    double* vec = (double *)malloc(vecSize * sizeof(double));
    vec[0] = 1.0;
    vec[1] = 2.0;
    vec[2] = 3.0;
    vec[3] = 4.0;
    vec[4] = 5.0;
    vec[5] = 6.0;

    NumericsMatrix * Arw_mat = Arrow_repr(vecSize, vec, varsCount);
    NumericsMatrix * D_Arw_mat = NM_create(NM_DENSE, vecSize, vecSize);
    NM_to_dense(Arw_mat, D_Arw_mat);
    NM_display(D_Arw_mat);

    assert(NM_get_value(Arw_mat, 0, 0) == 1.0);
    assert(NM_get_value(Arw_mat, 0, 1) == 2.0);
    assert(NM_get_value(Arw_mat, 0, 2) == 3.0);
    assert(NM_get_value(Arw_mat, 1, 0) == 2.0);
    assert(NM_get_value(Arw_mat, 2, 0) == 3.0);
    assert(NM_get_value(Arw_mat, 1, 1) == 1.0);
    assert(NM_get_value(Arw_mat, 2, 2) == 1.0);
    assert(NM_get_value(Arw_mat, 3, 3) == 4.0);
    assert(NM_get_value(Arw_mat, 3, 4) == 5.0);
    assert(NM_get_value(Arw_mat, 3, 5) == 6.0);
    assert(NM_get_value(Arw_mat, 4, 3) == 5.0);
    assert(NM_get_value(Arw_mat, 5, 3) == 6.0);
    assert(NM_get_value(Arw_mat, 4, 4) == 4.0);
    assert(NM_get_value(Arw_mat, 5, 5) == 4.0);
    return 0;
}

static int Arrow_repr_2d_test()
{
    int varsCount = 3;
    int dimension = 2;
    int vecSize = varsCount * dimension;

    double* vec = (double *)malloc(vecSize * sizeof(double));
    vec[0] = 1.0;
    vec[1] = 2.0;
    vec[2] = 3.0;
    vec[3] = 4.0;
    vec[4] = 5.0;
    vec[5] = 6.0;

    NumericsMatrix * Arw_mat = Arrow_repr(vecSize, vec, varsCount);
    NumericsMatrix * D_Arw_mat = NM_create(NM_DENSE, vecSize, vecSize);
    NM_to_dense(Arw_mat, D_Arw_mat);
    NM_display(D_Arw_mat);

    assert(NM_get_value(Arw_mat, 0, 0) == 1.0);
    assert(NM_get_value(Arw_mat, 0, 1) == 2.0);
    assert(NM_get_value(Arw_mat, 1, 0) == 2.0);
    assert(NM_get_value(Arw_mat, 1, 1) == 1.0);
    assert(NM_get_value(Arw_mat, 2, 2) == 3.0);
    assert(NM_get_value(Arw_mat, 2, 3) == 4.0);
    assert(NM_get_value(Arw_mat, 3, 2) == 4.0);
    assert(NM_get_value(Arw_mat, 3, 3) == 3.0);
    assert(NM_get_value(Arw_mat, 4, 4) == 5.0);
    assert(NM_get_value(Arw_mat, 4, 5) == 6.0);
    assert(NM_get_value(Arw_mat, 5, 4) == 6.0);
    assert(NM_get_value(Arw_mat, 5, 5) == 5.0);

    return 0;
}


int main(void)
{
    Arrow_repr_2d_test();
    Arrow_repr_3d_test();
    return 0;
}
