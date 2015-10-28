#!/usr/bin/env python

# Siconos-sample, Copyright INRIA 2005-2012.
# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.
# Siconos is a free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# Siconos is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Siconos; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
#
# Contact: Vincent ACARY, siconos-team@lists.gforge.fr

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
    except:
        pass
    if err_inv < 1.0e-12:
        try:
            (AexpS, PsiS) = compute_dt_matrices(A, B, h, TV)
            err_phi = np.linalg.norm((Aexp - AexpS))
            err_psi = np.linalg.norm(Psi - PsiS)
            if err_phi > 5.0e-12 or err_psi > 5.0e-12:
                print(err_phi, err_psi)
                failed = True
        except:
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
