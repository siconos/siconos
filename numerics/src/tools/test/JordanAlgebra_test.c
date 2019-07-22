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
#include "NumericsVector.h"
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


static int JA_eigenvals_test()
{
    const double EPS = 1e-8;
    int info = 0;
    int test_failed = 0;
    double out1[4];
    double vec[] = {1.256, 0.356, 0.874, 3.654, 0.154, 1.035};
    JA_eigenvals(&vec, 6, 2, &out1);

    test_failed += (fabs(out1[0] - 2.19972242) > EPS);
    test_failed += (fabs(out1[1] - 0.31227758) > EPS);
    test_failed += (fabs(out1[2] - 4.70039429) > EPS);
    test_failed += (fabs(out1[3] - 2.60760571) > EPS);

    if (test_failed > 0)
        info += 1;

    double out2[6];
    JA_eigenvals(&vec, 6, 3, &out2);

    test_failed += (fabs(out2[0] - 1.612) > EPS);
    test_failed += (fabs(out2[1] - 0.9) > EPS);
    test_failed += (fabs(out2[2] - 4.528) > EPS);
    test_failed += (fabs(out2[3] - (-2.78)) > EPS);
    test_failed += (fabs(out2[4] - 1.189) > EPS);
    test_failed += (fabs(out2[5] - (-0.881)) > EPS);


    if (test_failed > 0)
        info += 1;

    printf("== End of test JA_eigenvals_test(result = %d)\n", info);
    return info;
}


static int JA_eigenvecs_test()
{
    const double EPS = 1e-8;
    int info = 0;
    int test_failed = 0;
    double ** out1 = (double**)malloc(4 * sizeof(double*));
    for (int i = 0; i < 4; ++i)
        out1[i] = (double*)calloc(3, sizeof(double));

    double vec[] = {1.256, 0.356, 0.874, 3.654, 0.154, 1.035};
    JA_eigenvecs(&vec, 6, 2, out1);

    test_failed += (fabs(out1[0][0] - 0.5) > EPS);
    test_failed += (fabs(out1[0][1] - 0.18861478) > EPS);
    test_failed += (fabs(out1[0][2] - 0.46305989) > EPS);
    test_failed += (fabs(out1[1][0] - 0.5) > EPS);
    test_failed += (fabs(out1[1][1] + 0.18861478) > EPS);
    test_failed += (fabs(out1[1][2] + 0.46305989) > EPS);

    test_failed += (fabs(out1[2][0] - 0.5) > EPS);
    test_failed += (fabs(out1[2][1] - 0.07358603) > EPS);
    test_failed += (fabs(out1[2][2] - 0.49455545) > EPS);
    test_failed += (fabs(out1[3][0] - 0.5) > EPS);
    test_failed += (fabs(out1[3][1] + 0.07358603) > EPS);
    test_failed += (fabs(out1[3][2] + 0.49455545) > EPS);

    if (test_failed > 0)
        info += 1;

    for (int i = 0; i < 4; ++i)
        free(out1[i]);
    free(out1);

    out1 = (double**)malloc(6 * sizeof(double*));
    for (int i = 0; i < 6; ++i)
        out1[i] = (double*)calloc(2, sizeof(double));

    JA_eigenvecs(&vec, 6, 3, out1);

    test_failed += (fabs(out1[0][0] - 0.5) > EPS);
    test_failed += (fabs(out1[0][1] - 0.5) > EPS);
    test_failed += (fabs(out1[1][0] - 0.5) > EPS);
    test_failed += (fabs(out1[1][1] + 0.5) > EPS);
    test_failed += (fabs(out1[2][0] - 0.5) > EPS);
    test_failed += (fabs(out1[2][1] - 0.5) > EPS);
    test_failed += (fabs(out1[3][0] - 0.5) > EPS);
    test_failed += (fabs(out1[3][1] + 0.5) > EPS);
    test_failed += (fabs(out1[4][0] - 0.5) > EPS);
    test_failed += (fabs(out1[4][1] - 0.5) > EPS);
    test_failed += (fabs(out1[5][0] - 0.5) > EPS);
    test_failed += (fabs(out1[5][1] + 0.5) > EPS);


    if (test_failed > 0)
        info += 1;

    for (int i = 0; i < 6; ++i)
        free(out1[i]);
    free(out1);

    printf("== End of test JA_eigenvecs_test(result = %d)\n", info);
    return info;
}


static int JA_sqrt_test()
{
    const double EPS = 1e-8;
    int test_failed = 0;
    int info = 0;
    double vec[] = {1.256, 0.356, 0.874, 3.654, 0.154, 1.035};
    double * out = (double*)malloc(6 * sizeof(double));
    JA_sqrt(&vec, 6, 2, out);

    test_failed += (fabs(out[0] - 1.02098207) > EPS);
    test_failed += (fabs(out[1] - 0.17434194) > EPS);
    test_failed += (fabs(out[2] - 0.42801927) > EPS);
    test_failed += (fabs(out[3] - 1.89142377) > EPS);
    test_failed += (fabs(out[4] - 0.04071007) > EPS);
    test_failed += (fabs(out[5] - 0.27360341) > EPS);

    if (test_failed > 0)
        info += 1;

    free(out);
    printf("== End of test JA_sqrt_test(result = %d)\n", info);
    return info;
}


static int JA_det_test()
{
    const double EPS = 1e-8;
    int test_failed = 0;
    int info = 0;
    double vec[] = {1.256, 0.356, 0.874, 3.654, 0.154, 1.035};
    double * out = (double*)malloc(2 * sizeof(double));
    JA_det(&vec, 6, 2, out);

    test_failed += (fabs(out[0] - 0.686924) > EPS);
    test_failed += (fabs(out[1] - 12.256775) > EPS);

    if (test_failed > 0)
        info += 1;

    free(out);

    printf("== End of test JA_det_test(result = %d)\n", info);
    return info;
}

int main(void)
{
    int info = 0;
    info += Arrow_repr_2d_test();
    info += Arrow_repr_3d_test();
    info += JA_prod_test();
    info += JA_eigenvals_test();
    info += JA_eigenvecs_test();
    info += JA_sqrt_test();
    info += JA_det_test();

    return info;
}
