import os, sys, math
from collections import defaultdict

def generate(cM, output_dir):
    N =  (2 ** cM)
    os.makedirs(output_dir, exist_ok=True)
    cellSize = 1.0 / N

    node_map = {}
    next_nid = 1
    def get_nid(x, y):
        nonlocal next_nid
        key = (round(x, 14), round(y, 14))
        if key not in node_map:
            node_map[key] = next_nid
            next_nid += 1
        return node_map[key]

    for i in range(N + 1):
        for j in range(N + 1):
            get_nid(i / N, j / N)

    node_xyz = [None] * (next_nid)
    for (x, y), nid in node_map.items():
        node_xyz[nid] = (x, y)

    nEdgesPerEl = 4
    quads = []
    for i in range(N):
        for j in range(N):
            v00 = get_nid(i / N, j / N)
            v10 = get_nid((i + 1) / N, j / N)
            v11 = get_nid((i + 1) / N, (j + 1) / N)
            v01 = get_nid(i / N, (j + 1) / N)
            quads.append((v00, v10, v11, v01))

    nK = len(quads)
    all_tris_global = []
    edge_cnt = defaultdict(int)
    per_el_data = []

    for kidx, (v00, v10, v11, v01) in enumerate(quads):
        verts = [v00, v10, v11, v01]
        tri1 = (v00, v10, v11)
        tri2 = (v00, v11, v01)
        all_tris_global.extend([tri1, tri2])
        for (a, b, c) in [tri1, tri2]:
            for e2 in [(a, b), (b, c), (c, a)]:
                key = (min(e2[0], e2[1]), max(e2[0], e2[1]))
                edge_cnt[key] += 1

        local_of = {}
        local_coords = []
        for nid in verts:
            if nid not in local_of:
                local_of[nid] = len(local_of) + 1
                local_coords.append(node_xyz[nid])

        tris_local = [(local_of[a], local_of[b], local_of[c])
                      for (a, b, c) in [tri1, tri2]]

        bnd_local = []
        for e in range(4):
            label = 100 + nEdgesPerEl * kidx + e
            a, b = verts[e], verts[(e + 1) % 4]
            x0, y0 = node_xyz[a]
            x1, y1 = node_xyz[b]
            edge_len = math.hypot(x1 - x0, y1 - y0)
            nSub = max(1, int(edge_len / cellSize + 0.5))
            prev_lid = local_of[a]
            for s in range(1, nSub + 1):
                t = s / nSub
                xs = x0 + t * (x1 - x0)
                ys = y0 + t * (y1 - y0)
                nsid = get_nid(xs, ys)
                if nsid not in local_of:
                    local_of[nsid] = len(local_of) + 1
                    local_coords.append((xs, ys))
                new_lid = local_of[nsid]
                bnd_local.append((prev_lid, new_lid, label))
                prev_lid = new_lid

        per_el_data.append((len(local_coords), tris_local, bnd_local, local_coords))

    bnd_global = [e for e, cnt in edge_cnt.items() if cnt == 1]
    nnodes = len(node_map)
    node_xyz2 = [None] * (nnodes + 1)
    for (x, y), nid in node_map.items():
        node_xyz2[nid] = (x, y)

    eps = 1e-9
    with open(f"{output_dir}/calP.msh", "w") as f:
        f.write(f"{nnodes} {len(all_tris_global)} {len(bnd_global)}\n")
        for nid in range(1, nnodes + 1):
            x, y = node_xyz2[nid]
            f.write(f"{x:.15e} {y:.15e} 0\n")
        for (a, b, c) in all_tris_global:
            f.write(f"{a} {b} {c} 0\n")
        for (a, b) in bnd_global:
            xa, ya = node_xyz2[a]
            xb, yb = node_xyz2[b]
            xm, ym = (xa + xb) / 2, (ya + yb) / 2
            if abs(ym) < eps: lbl = 1
            elif abs(xm - 1) < eps: lbl = 2
            elif abs(ym - 1) < eps: lbl = 3
            elif abs(xm) < eps: lbl = 4
            else: lbl = 1
            f.write(f"{a} {b} {lbl}\n")

    print(f"calP: {nnodes} nodes, {len(all_tris_global)} tris, "
          f"{nK} quad elements, {len(bnd_global)} bnd edges")

    for kidx in range(nK):
        nv, tris, bnd, coords = per_el_data[kidx]
        with open(f"{output_dir}/Th_{kidx}.msh", "w") as f:
            f.write(f"{nv} {len(tris)} {len(bnd)}\n")
            for lid in range(nv):
                x, y = coords[lid]
                f.write(f"{x:.15e} {y:.15e} 0\n")
            for (a, b, c) in tris:
                f.write(f"{a} {b} {c} 0\n")
            for (a, b, lbl) in bnd:
                f.write(f"{a} {b} {lbl}\n")

    print(f"Generated {nK} Th_K.msh files in '{output_dir}/'")

if __name__ == "__main__":
    cM = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    out = f"meshes_quad_cM{cM}"
    generate(cM, out)
