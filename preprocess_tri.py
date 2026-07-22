import os, sys
from collections import Counter

# =====================================================================
# Triangular Macro Mesh Generator for the Full Unit Square [0,1]^2.
#
# Grid: 4*2^cM cells per side (same cell size as L-shape generator).
# Each cell is split into 2 triangles; each triangle is one macro
# element with 3 boundary edges.
#
# Level cM=0:  grid 4x4  → 32 elements,  h = 0.25
# Level cM=1:  grid 8x8  → 128 elements, h = 0.125
# Level cM=2:  grid 16x16 → 512 elements, h = 0.0625
#
# Labels: 100 + 3*K + e  (e = 0,1,2)
# Format: FreeFem++ native (.msh, readmesh)
# =====================================================================

def generate(cM, output_dir):
    N =  (2 ** cM)         # cells per side (matches L-shape cell size)
    h = 1.0 / N
    nn = N + 1
    os.makedirs(output_dir, exist_ok=True)

    # ---- Build global vertex grid ------------------------------------
    node_xyz = {}
    node_id = {}
    nid = 0
    for j in range(nn):
        for i in range(nn):
            nid += 1
            node_xyz[nid] = (i * h, j * h)
            node_id[(i, j)] = nid
    nnodes = nid

    def gid(i, j):
        return node_id[(i, j)]

    def vertex_label(i, j):
        """Domain boundary label for vertex (i,j) in [0,nn-1]^2."""
        if i == 0 and j == 0:        return 1
        if i == 0 and j == nn - 1:   return 4
        if i == nn - 1 and j == 0:   return 1
        if i == nn - 1 and j == nn - 1: return 3
        if j == 0:                   return 1
        if i == nn - 1:              return 2
        if j == nn - 1:              return 3
        if i == 0:                   return 4
        return 0

    nEdgesPerEl = 3
    total_elements = 2 * N * N

    elem_data = []  # [(K, vlist_global, tris_local, bnd_local)]

    for j in range(N):
        for i in range(N):
            v00, v10, v11, v01 = (i, j), (i + 1, j), (i + 1, j + 1), (i, j + 1)

            # ---- Element K0 = (i,j)-(i+1,j)-(i+1,j+1) ----
            K0 = 2 * (j * N + i)
            # CCW boundary order
            bnd_ordered0 = [v00, v10, v11]
            vlist0 = bnd_ordered0  # no interior vertices
            local0 = {v: idx + 1 for idx, v in enumerate(vlist0)}
            tris_local0 = [tuple(local0[v] for v in [v00, v10, v11])]
            bnd_local0 = []
            for e in range(3):
                va = bnd_ordered0[e]
                vb = bnd_ordered0[(e + 1) % 3]
                label = 100 + nEdgesPerEl * K0 + e
                bnd_local0.append((local0[va], local0[vb], label))
            elem_data.append((K0, vlist0, tris_local0, bnd_local0))

            # ---- Element K1 = (i,j)-(i+1,j+1)-(i,j+1) ----
            K1 = 2 * (j * N + i) + 1
            bnd_ordered1 = [v00, v11, v01]
            vlist1 = bnd_ordered1
            local1 = {v: idx + 1 for idx, v in enumerate(vlist1)}
            tris_local1 = [tuple(local1[v] for v in [v00, v11, v01])]
            bnd_local1 = []
            for e in range(3):
                va = bnd_ordered1[e]
                vb = bnd_ordered1[(e + 1) % 3]
                label = 100 + nEdgesPerEl * K1 + e
                bnd_local1.append((local1[va], local1[vb], label))
            elem_data.append((K1, vlist1, tris_local1, bnd_local1))

    # ---- Global triangulation for calP.msh ---------------------------
    all_tris_global = []
    for j in range(N):
        for i in range(N):
            a, b, c, d = gid(i, j), gid(i + 1, j), gid(i + 1, j + 1), gid(i, j + 1)
            all_tris_global.append((a, b, c))
            all_tris_global.append((a, c, d))

    # Boundary edges: edges appearing only once
    edge_cnt = Counter()
    for (a, b, c) in all_tris_global:
        for e in [(a, b), (b, c), (c, a)]:
            key = e if e[0] < e[1] else (e[1], e[0])
            edge_cnt[key] += 1
    bnd_global_edges = [e for e, cnt in edge_cnt.items() if cnt == 1]

    # Sort boundary edges for FreeFem++ extract() compatibility
    def bnd_sort_key(item):
        a, b = item
        xa, ya = node_xyz[a]
        xb, yb = node_xyz[b]
        xm, ym = (xa + xb) / 2, (ya + yb) / 2
        if ya == 0.0 and yb == 0.0:
            return (1, xm)
        elif xa == 1.0 and xb == 1.0:
            return (2, -ym)
        elif ya == 1.0 and yb == 1.0:
            return (3, -xm)
        elif xa == 0.0 and xb == 0.0:
            return (4, ym)
        else:
            return (0, xm)

    bnd_global_edges_sorted = sorted(bnd_global_edges, key=bnd_sort_key)

    # ---- Write calP.msh (FreeFem++ native format) --------------------
    with open(f"{output_dir}/calP.msh", "w") as f:
        f.write(f"{nnodes} {len(all_tris_global)} {len(bnd_global_edges)}\n")
        for n in sorted(node_xyz.keys()):
            x, y = node_xyz[n]
            lbl = vertex_label(*[(i, j) for i, j in node_id if node_id[(i, j)] == n][0])
            f.write(f"{x:.15e} {y:.15e} {lbl}\n")
        for (a, b, c) in all_tris_global:
            f.write(f"{a} {b} {c} 0\n")
        for (a, b) in bnd_global_edges_sorted:
            xa, ya = node_xyz[a]
            xb, yb = node_xyz[b]
            if ya == 0.0 and yb == 0.0:
                lbl = 1
            elif ya == 1.0 and yb == 1.0:
                lbl = 3
            elif xa == 0.0 and xb == 0.0:
                lbl = 4
            elif xa == 1.0 and xb == 1.0:
                lbl = 2
            else:
                lbl = 1
            f.write(f"{a} {b} {lbl}\n")

    print(f"calP: {nnodes} nodes, {len(all_tris_global)} tris, "
          f"{total_elements} tri elements, {len(bnd_global_edges)} bnd edges")

    # ---- Write Th_K.msh (FreeFem++ native format) --------------------
    for (K, vlist, tris_local, bnd_local) in elem_data:
        with open(f"{output_dir}/Th_{K}.msh", "w") as f:
            nv = len(vlist)
            ntri = len(tris_local)
            nbnd = len(bnd_local)
            f.write(f"{nv} {ntri} {nbnd}\n")
            for lid, (gi, gj) in enumerate(vlist, 1):
                f.write(f"{gi * h:.15e} {gj * h:.15e} 0\n")
            for (a, b, c) in tris_local:
                f.write(f"{a} {b} {c} 0\n")
            for (a, b, lbl) in bnd_local:
                f.write(f"{a} {b} {lbl}\n")

    print(f"Generated {len(elem_data)} Th_K.msh files in '{output_dir}/'")


if __name__ == "__main__":
    cM = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    out = f"meshes_tri_cM{cM}"
    generate(cM, out)
