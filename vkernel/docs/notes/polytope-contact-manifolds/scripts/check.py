"""Compare, step by step, the active set predicted from the manifold port
(manifold.py) and the Moreau-Jean activation rule with the simulator's.

   python3 check.py <scene>..."""

import numpy as np, sys
from sim import history

for s in sys.argv[1:]:
    ts, q, v, cf, hist, worst, mismatch = history(s)
    print(
        f"{s}: steps where predicted active-set size != simulator's: {mismatch}/{len(ts) - 1}; cp match worst {worst:.1e}"
    )
    for label, hh in hist.items():
        a = hh["active"]
        print(
            f"   {label}: active steps p0 {a[:, 0].sum()} p1 {a[:, 1].sum()} both {(a[:, 0] & a[:, 1]).sum()}"
        )
