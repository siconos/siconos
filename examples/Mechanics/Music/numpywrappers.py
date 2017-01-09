"""Interface to numpy arrays

Those functions are useful to enforce predefined types in numpy arrays
and avoid unexpected behavior due to types conversion.
"""
import numpy as np


REAL_KIND = np.float64
"""Set default type for real numbers"""

ORDER = 'F'
"""default array layout (fortran or C convention)"""


def zeros(shape, dtype=REAL_KIND):
    """
    Wrapper to numpy.zeros, force order to hysop.constants.ORDER
    """
    return np.zeros(shape, dtype=dtype, order=ORDER)


def ones(shape, dtype=REAL_KIND):
    """
    Wrapper to numpy.ones, force order to hysop.constants.ORDER
    """
    return np.ones(shape, dtype=dtype, order=ORDER)


def zeros_like(tab):
    """
    Wrapper to numpy.zeros_like, force order to hysop.constants.ORDER
    """
    return np.zeros_like(tab, dtype=tab.dtype, order=ORDER)


def asarray(tab):
    """
    Wrapper to numpy.asarray, force order to hysop.constants.ORDER
    """
    return np.asarray(tab, order=ORDER, dtype=tab.dtype)


def asrealarray(tab):
    """
    Wrapper to numpy.asarray, force order to hysop.constants.ORDER
    and type to hysop.constants.HYSOP_REAL
    """
    return np.asarray(tab, order=ORDER, dtype=REAL_KIND)


# Some useful functions, from http://ipython-books.github.io/featured-01/
def get_data_base(arr):
    """For a given Numpy array, finds the
    base array that "owns" the actual data."""
    base = arr
    while isinstance(base.base, np.ndarray):
        base = base.base
    return base


def arrays_share_data(x, y):
    """Returns true if x and y arrays share the same
    memory location
    """
    return get_data_base(x) is get_data_base(y)
