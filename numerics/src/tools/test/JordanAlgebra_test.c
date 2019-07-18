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
#include "math.h"

static int Arrow_repr_3d_test()
{
    int info = 0;
    int test_failed = 0;
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

    NumericsMatrix * Arw_mat = Arrow_repr(vec, vecSize, varsCount);
    NumericsMatrix * D_Arw_mat = NM_create(NM_DENSE, vecSize, vecSize);
    NM_to_dense(Arw_mat, D_Arw_mat);
    // NM_display(D_Arw_mat);

    test_failed += !(NM_get_value(Arw_mat, 0, 0) == 1.0);
    test_failed += !(NM_get_value(Arw_mat, 0, 1) == 2.0);
    test_failed += !(NM_get_value(Arw_mat, 0, 2) == 3.0);
    test_failed += !(NM_get_value(Arw_mat, 1, 0) == 2.0);
    test_failed += !(NM_get_value(Arw_mat, 2, 0) == 3.0);
    test_failed += !(NM_get_value(Arw_mat, 1, 1) == 1.0);
    test_failed += !(NM_get_value(Arw_mat, 2, 2) == 1.0);
    test_failed += !(NM_get_value(Arw_mat, 3, 3) == 4.0);
    test_failed += !(NM_get_value(Arw_mat, 3, 4) == 5.0);
    test_failed += !(NM_get_value(Arw_mat, 3, 5) == 6.0);
    test_failed += !(NM_get_value(Arw_mat, 4, 3) == 5.0);
    test_failed += !(NM_get_value(Arw_mat, 5, 3) == 6.0);
    test_failed += !(NM_get_value(Arw_mat, 4, 4) == 4.0);
    test_failed += !(NM_get_value(Arw_mat, 5, 5) == 4.0);
    if (test_failed > 0)
        info = 1;

    printf("== End of test Arrow_repr_3d_test(result = %d)\n", info);
    return info;
}

static int Arrow_repr_2d_test()
{
    int info = 0;
    int test_failed = 0;
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

    NumericsMatrix * Arw_mat = Arrow_repr(vec, vecSize, varsCount);
    NumericsMatrix * D_Arw_mat = NM_create(NM_DENSE, vecSize, vecSize);
    NM_to_dense(Arw_mat, D_Arw_mat);
    // NM_display(D_Arw_mat);

    test_failed += !(NM_get_value(Arw_mat, 0, 0) == 1.0);
    test_failed += !(NM_get_value(Arw_mat, 0, 1) == 2.0);
    test_failed += !(NM_get_value(Arw_mat, 1, 0) == 2.0);
    test_failed += !(NM_get_value(Arw_mat, 1, 1) == 1.0);
    test_failed += !(NM_get_value(Arw_mat, 2, 2) == 3.0);
    test_failed += !(NM_get_value(Arw_mat, 2, 3) == 4.0);
    test_failed += !(NM_get_value(Arw_mat, 3, 2) == 4.0);
    test_failed += !(NM_get_value(Arw_mat, 3, 3) == 3.0);
    test_failed += !(NM_get_value(Arw_mat, 4, 4) == 5.0);
    test_failed += !(NM_get_value(Arw_mat, 4, 5) == 6.0);
    test_failed += !(NM_get_value(Arw_mat, 5, 4) == 6.0);
    test_failed += !(NM_get_value(Arw_mat, 5, 5) == 5.0);

    if (test_failed > 0)
        info = 1;

    printf("== End of test Arrow_repr_2d_test(result = %d)\n", info);
    return info;
}

static int JA_prod_test()
{
    int info = 0;
    int test_failed = 0;
    int n = 6;
    const double EPS = 1e-12;
    double x[] = {0.75, 0.35, 0.1, 0.78, 1.25, 0.78};
    double y[] = {0.27, 4.35, 1.02, 0.35, 0.78, 0.236};

    double * res = (double*)calloc(n, sizeof(double));
    JA_prod(&x, &y, n, 1, res);
    test_failed += (fabs(res[0] - 3.25908) > EPS);
    test_failed += (fabs(res[1] - 3.357) > EPS);
    test_failed += (fabs(res[2] - 0.792) > EPS);
    test_failed += (fabs(res[3] - 0.4731) > EPS);
    test_failed += (fabs(res[4] - 0.9225) > EPS);
    test_failed += (fabs(res[5] - 0.3876) > EPS);

    if (test_failed > 0)
        info += 1;

    test_failed = 0;

    JA_prod(&x, &y, n, 2, res);
    test_failed += (fabs(res[0] - 1.827) > EPS);
    test_failed += (fabs(res[1] - 3.357) > EPS);
    test_failed += (fabs(res[2] - 0.792) > EPS);
    test_failed += (fabs(res[3] - 1.43208) > EPS);
    test_failed += (fabs(res[4] - 1.0459) > EPS);
    test_failed += (fabs(res[5] - 0.45708) > EPS);

    if (test_failed > 0)
        info += 1;

    printf("== End of test JA_prod_test(result = %d)\n", info);
    free(res);
    return info;
}


int main(void)
{
    int info = 0;
    info += Arrow_repr_2d_test();
    info += Arrow_repr_3d_test();
    info += JA_prod_test();
    return info;
}
