# The unit square as a gmsh 2.2 mesh of N x N cells, two T3 per cell, with
# the physical groups of mesh_N5.msh (N = 5): the bulk material (1), the
# Dirichlet nodes at the bottom (2), an applied force node at the top left
# (3), the contact nodes all around (4, segments counter-clockwise from the
# bottom left node).
#
# usage: python3 square_mesh.py <N> <output.msh>

import sys

n = int(sys.argv[1])
output = sys.argv[2]


def node(i, j):
    """the node of column i, row j (from 1)"""
    return j * (n + 1) + i + 1


elements = []
for j in range(n):
    for i in range(n):
        a, b = node(i, j), node(i + 1, j)
        c, d = node(i + 1, j + 1), node(i, j + 1)
        elements.append(f"2 2 1 7 {a} {b} {c}")
        elements.append(f"2 2 1 7 {a} {c} {d}")

# the boundary, counter-clockwise from the bottom left node
boundary = (
    [node(i, 0) for i in range(n + 1)]
    + [node(n, j) for j in range(1, n + 1)]
    + [node(i, n) for i in range(n - 1, -1, -1)]
    + [node(0, j) for j in range(n - 1, -1, -1)]
)
for a, b in zip(boundary[:-1], boundary[1:]):
    elements.append(f"1 2 4 4 {a} {b}")

for i in range(n + 1):
    elements.append(f"15 2 2 2 {node(i, 0)}")
elements.append(f"15 2 3 3 {node(0, n)}")

with open(output, "w") as f:
    f.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n")
    f.write("$PhysicalNames\n4\n")
    f.write('2 1 "Bulk Material"\n1 2 "Dirichlet BC"\n0 3 "Applied Force"\n')
    f.write('1 4 "Contact Nodes"\n$EndPhysicalNames\n')
    f.write(f"$Nodes\n{(n + 1) ** 2}\n")
    for j in range(n + 1):
        for i in range(n + 1):
            f.write(f"{node(i, j)} {i / n:.16e} {j / n:.16e} {0.0:.16e}\n")
    f.write("$EndNodes\n")
    f.write(f"$Elements\n{len(elements)}\n")
    for k, e in enumerate(elements, start=1):
        f.write(f"{k} {e}\n")
    f.write("$EndElements\n")
