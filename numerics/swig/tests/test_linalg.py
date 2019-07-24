import siconos.numerics as sn
import scipy.sparse
import numpy as np
import copy
from siconos.tests_setup import working_dir


data_dir = working_dir + '/data/'

dbl_eps = np.finfo(np.float).eps


def check_size(sa, sb):
    assert sa[0] == sb[0]
    assert sa[1] == sb[1]


def compare_with_SBM(sbmat, mat):
    print('in compare_with_SBM')
    if scipy.sparse.issparse(mat):
        mdense = mat.todense()

    m, n = mdense.shape

    for i in range(m):
        for j in range(n):
            assert sn.SBM_get_value(sbmat, i, j) == mdense[i, j]


def SBM_tests(sbm):
    print('in SBM_tests')
    # SBM_to_dense
    # SBM_to_sparse
    # SBM_from_csparse
    nm = sn.NumericsMatrix(sbm)
    assert nm is not None
    mcoo = sn.NM_triplet(nm)
    mcoo._check()
    check_size((nm.size0, nm.size1), mcoo.shape)
    compare_with_SBM(sbm, mcoo)

    mcsc = sn.NM_csc(nm)
    mcsc.check_format(True)
    check_size((nm.size0, nm.size1), mcsc.shape)
    compare_with_SBM(sbm, mcsc)

    mcsc_trans = sn.NM_csc_trans(nm)
    mcsc_trans.check_format(True)
    check_size((nm.size0, nm.size1), mcsc_trans.T.shape)
    compare_with_SBM(sbm, mcsc_trans.T)

    return copy.deepcopy(mcsc)


def dense_tests(mdense):
    print('in dense_tests')
    nm = sn.NumericsMatrix(mdense)
    check_size((nm.size0, nm.size1), mdense.shape)
    assert nm is not None


def sparse_tests(spm):
    print('in sparse_tests')
#    sp_formats = ('bsr', 'coo', 'csc', 'csr', 'dia', 'dok', 'lil')
    sp_formats = ('bsr', 'coo', 'csc', 'csr', 'dia', 'lil')
    spm.check_format(True)

    spm_shape = spm.shape
    nm = sn.NumericsMatrix(spm)
    assert nm is not None
    check_size((nm.size0, nm.size1), spm_shape)

    for fmt in sp_formats:
        func = getattr(scipy.sparse, fmt + '_matrix')
        mm = func(spm)
        if fmt in ('csc', 'csr'):
            mm.check_format(True)
        elif fmt == 'coo':
            mm._check()

        nm = sn.NumericsMatrix(mm)
        check_size((nm.size0, nm.size1), spm_shape)
        assert nm is not None
        t1 = sn.NM_triplet(nm)
        t1._check()
        check_size(t1.shape, spm_shape)
        assert ((spm - t1).nnz == 0)

        t2 = sn.NM_csc(nm)
        t2.check_format(True)
        check_size(t2.shape, spm_shape)
        assert ((spm - t2).nnz == 0)

        t2bis = sn.NM_csc_trans(nm)
        t2bis.check_format(True)
        check_size(t2bis.T.shape, spm_shape)
        assert ((spm.T - t2bis).nnz == 0)

        t3 = sn.NM_triplet(mm)
        t3._check()
        check_size(t3.shape, spm_shape)
        assert ((spm - t3).nnz == 0)

        t4 = sn.NM_csc(mm)
        t4.check_format(True)
        check_size(t4.shape, spm_shape)
        assert ((spm - t4).nnz == 0)

        t4bis = sn.NM_csc_trans(mm)
        t4bis.check_format(True)
        check_size(t4bis.T.shape, spm_shape)
        assert ((spm.T - t4bis).nnz == 0)


def test_create():
    fcp = sn.GlobalFrictionContactProblem()

    Min = np.ones((12, 3))
    Hin = np.ones((3, 3))

    fcp.M = Min
    fcp.H = Hin

    assert np.max(Min - fcp.M) == 0.
    assert np.max(Hin - fcp.H) == 0.

    sp_formats = ('bsr', 'coo', 'csc', 'csr', 'dia', 'lil')

    for fmt in sp_formats:
        func = getattr(scipy.sparse, fmt + '_matrix')
        mm = func(Min)
        hh = func(Hin)

        fcp.M = mm
        fcp.H = hh

        assert np.max(Min - fcp.M.todense()) == 0.
        assert np.max(Hin - fcp.H.todense()) == 0.


def test_convert():

    mat = []

    fcp = sn.GlobalFrictionContactProblem()
    sn.globalFrictionContact_newFromFile(fcp, data_dir + 'GFC3D_TwoRods1.dat')

    mat.append(fcp.M)
    mat.append(fcp.H)

    hacklist = [fcp]

    try:
        # cheap test ...
        data_fclib = ('LMGC_GFC3D_CubeH8.hdf5', 'LMGC_GlobalFrictionContactProblem00046.hdf5')

        for d in data_fclib:
            fcp = sn.globalFrictionContact_fclib_read(data_dir + d)
            hacklist.append(fcp)

            # Ho LMGC
            Mdense = fcp.M.todense()
            Mdense[np.nonzero(Mdense <= 100*dbl_eps)] = 0.
            MM = scipy.sparse.coo_matrix(Mdense)

            Hdense = fcp.H.todense()
            Hdense[np.nonzero(Hdense <= 100*dbl_eps)] = 0.
            HH = scipy.sparse.coo_matrix(Hdense)

            mat.append(MM)
            mat.append(HH)
    except:
        pass

    for m in mat:
        print("testing")
        if not scipy.sparse.issparse(m) and not isinstance(m, np.ndarray):
            spm = SBM_tests(m)
        else:
            spm = scipy.sparse.csc_matrix(m)

        mdense = spm.todense()
        dense_tests(mdense)
        sparse_tests(spm)


def test_arith():
    # NM_gemv
    pass


def test_solve():
    # NM_gesv
    pass
