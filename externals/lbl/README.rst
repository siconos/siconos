===================================
LBL: C wrapper around MA27 and MA57
===================================

:Author: Dominique Orban <dominique.orban@gerad.ca>

LBL is a unified C interface to the multifrontal symmetric indefinite linear
system solvers MA27 and MA57 from the Harwell Subroutine Library.

Example
=======

Here are the essential bits of a typical example (see `examples/examples.c`
for the complete code)::

    #include "lbl.h"
    LBL_Data *lbl;           // Main data structure.
    lbl = LBL_Initialize(nnz, n, stderr, 1);
    LBL_Analyze(lbl, 0);     // 0 = automatic pivot choice.
    LBL_Factorize(lbl, val); // val is an array of doubles.
    LBL_Solve(lbl, rhs);     // rhs is an array of doubles.
    LBL_Finalize(lbl);

Installing LBL
==============

1. Place the source files for MA27 and MA57 in some directory. You will need
   fd15d.f ma27ad.f mc21d.f mc34d.f mc47d.f mc59d.f mc64d.f mc71d.f ma57d.f.
   See http://www.hsl.rl.ac.uk/hsl2007

2. Edit and customize make.inc

3. Type ``make all``.

You can clean out temporary files with ``make clean`` and sanitize everything
with ``make purge``.


Testing LBL
===========

Type ``make example``. Change to the `example` directory and type
``./example``.

The expected output is::

   Input matrix:
   1  1  2.000000
   1  2  3.000000
   2  3  4.000000
   2  5  6.000000
   3  3  1.000000
   3  4  5.000000
   5  5  1.000000
   Input right-hand side:
   8.000000 45.000000 31.000000 15.000000 17.000000 
   Solution:
   1.000000 2.000000 3.000000 4.000000 5.000000


Bugs, Comments, Feature Requests
================================

Please send bugs, comments, feature requests, test cases, chocolate and coffee
to the `Lighthouse LBL page
<http://pykrylov.lighthouseapp.com/projects/54633-lbl>`_.


Enjoy!


.. image:: https://d2weczhvl823v0.cloudfront.net/dpo/lbl/trend.png
   :alt: Bitdeli badge
   :target: https://bitdeli.com/free

