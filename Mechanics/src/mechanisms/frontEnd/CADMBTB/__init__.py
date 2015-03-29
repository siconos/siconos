from Siconos.cadmbtb import cadmbtb
from Siconos.cadmbtb.cadmbtb import *

_imported_modules = [cadmbtb]

_bad_items = ['SHARED_PTR_DISOWN', 'os', 'sys']
_bad_expr = ['_swigregister', 'nullDeleter', 'weakref']

__all__ = []

for mod in _imported_modules:
    _list_items = dir(mod)

    for _it in _list_items:
        _not_add = 0
        if _it in _bad_items or _it[0] == '_':
            break

        for _bad in _bad_expr:
            if _bad in _it:
                _not_add = 1
                break

        if _not_add == 0:
            __all__.append(_it)
