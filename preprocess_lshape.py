import os, sys

# =====================================================================
# L-Shaped Macro Mesh Generator for the Full Unit Square [0,1]^2.
# 
# Each L-element contains exactly 6 corner vertices (no intermediate nodes)
# and is split into EXACTLY 4 triangles radiating from its re-entrant corner.
#
# FreeFem++ BAMG Safe: Enforces strict CCW orientation and contiguous 
# global boundary loops to prevent non-manifold topology errors.
# =====================================================================

# Ordered 6 corners for each L-shape type relative to the macro-block corner
# These follow a strict Counter-Clockwise (CCW) perimeter path.
CORNERS = [
    [(0,0), (2,0), (2,3), (1,3), (1,1), (0,1)], # gi=0
    [(0,1), (1,1), (1,3), (2,3), (2,4), (0,4)], # gi=1
    [(4,4), (2,4), (2,3), (3,3), (3,1), (4,1)], # gi=2
    [(4,1), (3,1), (3,3), (2,3), (2,0), (4,0)]  # gi=3
]

# 4-triangle local connectivity patterns (0-indexed based on CORNERS array)
# Triangles radiate outward from the unique re-entrant corner of each L-shape.
TRI_PATTERNS = [
    [(4, 5, 0), (4, 0, 1), (4, 1, 2), (4, 2, 3)], # gi=0 (radiates from index 4: (1,1))
    [(2, 3, 4), (2, 4, 5), (2, 5, 0), (2, 0, 1)], # gi=1 (radiates from index 2: (1,3))
    [(3, 4, 5), (3, 5, 0), (3, 0, 1), (3, 1, 2)], # gi=2 (radiates from index 3: (3,3))
    [(1, 2, 3), (1, 3, 4), (1, 4, 5), (1, 5, 0)]  # gi=3 (radiates from index 1: (3,1))
]

def vertex_label(i, j, nb):
    """Domain boundary label for vertex (i,j) in [0,nb]^2."""
    if i == 0 and j == 0:        return 1
    if i == 0 and j == nb:       return 4
    if i == nb and j == 0:       return 1
    if i == nb and j == nb:      return 3
    if j == 0:                   return 1
    if i == nb:                  return 2
    if j == nb:                  return 3
    if i == 0:                   return 4
    return 0

def generate(cM, output_dir):
    nb = 4 * (2**cM)
    h = 1.0 / nb
    os.makedirs(output_dir, exist_ok=True)

    nbx = nb // 4

    # ---- Step 1: Collect strictly unique global corner coordinates ----
    global_coords = set()
    for by in range(nbx):
        for bx in range(nbx):
            i0 = 4 * bx
            j0 = 4 * by
            for gi in range(4):
                for (cx, cy) in CORNERS[gi]:
                    global_coords.add((i0 + cx, j0 + cy))

    # Sort globally to create stable node numbering
    sorted_coords = sorted(list(global_coords), key=lambda p: (p[1], p[0]))
    node_id = {coords: idx for idx, coords in enumerate(sorted_coords, 1)}
    id_to_ij = {idx: coords for coords, idx in node_id.items()}
    nnodes = len(node_id)

    all_tris_global = []
    elem_data = []

    # ---- Step 2: Build Elements ----
    for by in range(nbx):
        for bx in range(nbx):
            i0 = 4 * bx
            j0 = 4 * by

            for gi in range(4):
                K = gi + 4 * (bx + nbx * by)
                corners = CORNERS[gi]
                tri_pattern = TRI_PATTERNS[gi]

                # Map local corners to global integer grid coordinates
                vlist = [(i0 + cx, j0 + cy) for (cx, cy) in corners]
                
                tris_local = []
                tris_global = []
                for tri in tri_pattern:
                    # Map to 1-based local indexing for FreeFem++
                    tris_local.append((tri[0]+1, tri[1]+1, tri[2]+1))
                    # Map to global IDs
                    tris_global.append((node_id[vlist[tri[0]]], node_id[vlist[tri[1]]], node_id[vlist[tri[2]]]))
                    
                all_tris_global.extend(tris_global)

                # Boundary sub-edges (exactly 6 un-cut edges per L-element)
                bnd_local = []
                for eidx in range(6):
                    label = 100 + 6*K + eidx
                    # Consecutive CCW local edges: (1,2), (2,3), (3,4), (4,5), (5,6), (6,1)
                    va_local = eidx + 1
                    vb_local = (eidx + 1) % 6 + 1
                    bnd_local.append((va_local, vb_local, label))

                elem_data.append((K, vlist, tris_local, bnd_local))

    # ---- Step 3: Global Boundary Construction ----
    # BAMG strictly requires a contiguous CCW closed boundary loop. 
    # To prevent any geometry-sorting errors, we manually walk the 4 walls of the domain.
    calP_bnd_edges = []
    
    # Bottom Edge (y=0), walking Left to Right
    bot = sorted([n for n, (i,j) in id_to_ij.items() if j == 0], key=lambda n: id_to_ij[n][0])
    for k in range(len(bot)-1): calP_bnd_edges.append((bot[k], bot[k+1], 1))
        
    # Right Edge (x=nb), walking Bottom to Top
    right = sorted([n for n, (i,j) in id_to_ij.items() if i == nb], key=lambda n: id_to_ij[n][1])
    for k in range(len(right)-1): calP_bnd_edges.append((right[k], right[k+1], 2))
        
    # Top Edge (y=nb), walking Right to Left
    top = sorted([n for n, (i,j) in id_to_ij.items() if j == nb], key=lambda n: id_to_ij[n][0], reverse=True)
    for k in range(len(top)-1): calP_bnd_edges.append((top[k], top[k+1], 3))
        
    # Left Edge (x=0), walking Top to Bottom
    left = sorted([n for n, (i,j) in id_to_ij.items() if i == 0], key=lambda n: id_to_ij[n][1], reverse=True)
    for k in range(len(left)-1): calP_bnd_edges.append((left[k], left[k+1], 4))

    # ---- Step 4: Write calP.msh ----
    with open(os.path.join(output_dir, "calP.msh"), "w") as f:
        f.write(f"{nnodes} {len(all_tris_global)} {len(calP_bnd_edges)}\n")
        for n in range(1, nnodes + 1):
            i, j = id_to_ij[n]
            x, y = i * h, j * h
            lbl = vertex_label(i, j, nb)
            f.write(f"{x:.15e} {y:.15e} {lbl}\n")
        for (a, b, c) in all_tris_global:
            f.write(f"{a} {b} {c} 0\n")
        for (a, b, lbl) in calP_bnd_edges:
            f.write(f"{a} {b} {lbl}\n")

    n_actual = len(elem_data)
    print(f"calP: {nnodes} nodes, {len(all_tris_global)} tris, {n_actual} L-elements, {len(calP_bnd_edges)} bnd edges")

    # ---- Step 5: Write Individual Submeshes Th_K.msh ----
    for (K, vlist, tris_local, bnd_local) in elem_data:
        with open(os.path.join(output_dir, f"Th_{K}.msh"), "w") as f:
            nv = 6
            ntri = 4
            nbnd = 6
            f.write(f"{nv} {ntri} {nbnd}\n")
            for (gi, gj) in vlist:
                x = gi * h
                y = gj * h
                f.write(f"{x:.15e} {y:.15e} 0\n")
            for (a, b, c) in tris_local:
                f.write(f"{a} {b} {c} 0\n")
            for (a, b, lbl) in bnd_local:
                f.write(f"{a} {b} {lbl}\n")

    print(f"Generated {n_actual} Th_K.msh files in '{output_dir}/'")

if __name__ == "__main__":
    cM = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    out = f"meshes_lshape_cM{cM}"
    generate(cM, out)