#!/usr/bin/env python

# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.
#
# Copyright 2016 INRIA.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import siconos.kernel as SK
import siconos.numerics as SN
import numpy as np
from scipy.linalg import expm
from siconos.functions import compute_dt_matrices

kmax = 10 if SN.compiled_in_debug_mode() else 100


def check_error(n, m, h, TV=False):
    err_inv = 1
    failed = False
    A = np.random.random((n, n))
    B = np.random.random((n, m))
    try:
        Ainv = np.linalg.inv(A)
        Aexp = expm(A * h)
        Psi = Ainv.dot((Aexp - np.eye(n)).dot(B))
        err_inv = np.linalg.norm(A.dot(np.linalg.inv(A)) - np.eye(n))
    except Exception as e :
        print('Exception in check_error(I) :', e)
        pass

    print 'err_inv:', err_inv
    if err_inv < 1.0e-12:
        try:
            (AexpS, PsiS) = compute_dt_matrices(A, B, h, TV)
            err_phi = np.linalg.norm((Aexp - AexpS))
            err_psi = np.linalg.norm(Psi - PsiS)
            if err_phi > 5.0e-12 or err_psi > 5.0e-12:
                print(err_phi, err_psi)
                failed = True
        except Exception as e :
            print('Exception in check_error(II) :', e)
            print(A.shape)
            failed = True
    return failed



def test_TI():
    h = 1.0e-4
    k = 0
    err = False
    while k < kmax:
        n = np.random.randint(2, kmax)
        m = np.random.randint(1, n)
        err = check_error(n, m, h)
        k += 1

    assert not err


def test_SISO():
    h = 1.0e-4
    k = 0
    err = False
    while k < kmax:
        n = np.random.randint(2, kmax)
        m = 1
        err = check_error(n, m, h)
        k += 1

    assert not err



def test_TV():
    h = 1.0e-4
    k = 0
    err = False
    while k < 10:
        n = np.random.randint(2, kmax)
        m = np.random.randint(1, n)
        err = check_error(n, m, h, True)
        k += 1

    assert not err

if __name__ == '__main__':
    print 'test_TI'
    test_TI()
