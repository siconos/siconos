from IO import *

__all__ = filter(lambda s:
                 any(b not in s for b in
                     ['_swigregister', 'nullDeleter', 'weakref']),
                 filter(lambda s:
                        s[0] != '_' and
                        s not in ['Kernel', 'SHARED_PTR_DISOWN', 'os', 'sys'],
                        dir(IO)))

