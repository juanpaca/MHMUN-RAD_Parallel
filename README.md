# MHM-RAD Parallel

MPI-parallel copy of the optimized polytopal MHM-RAD solver.

## Files

- `MHM-RAD_parallel.edp` - dispatcher for the supported polynomial pairs.
- `MHM-RAD_parallel_body.idp` - MPI version of the solver body.

## Run

```bash
ff-mpirun -np 4 MHM-RAD_parallel.edp -cM 1 -edgeM 0 -subM 0 -EkO 1 -PkO 3 -stab 1 -verboseinK 0
```

Default meshes are read from `../meshes_lshape_cM<N>/`. Override with `-meshDir <path>`. The meshes are generated with the python scripts.

## Parallelization

- Each rank builds the common mesh and skeleton metadata.
- W0 detection is distributed over contiguous macro-element blocks and reduced with MPI.
- Offline local problems are distributed over the same macro-element blocks.
- Local Schur triplets are gathered to rank 0.
- Offline reconstruction arrays are reduced so rank 0 can solve and postprocess exactly as the serial reference.

## Validation Against `MHM-RAD_poly_opt.edp`

Representative L-shaped mesh results matched the current optimized serial solver:

| Case | MPI ranks | Result |
|---|---:|---|
| `cM=0 edgeM=0 L1K3 stab=1` | 2 | `0.75 0.75 0.75 32 0.0282032 0.493064` |
| `cM=1 edgeM=0 L1K3 stab=1` | 2 | `0.375 0.375 0.375 112 0.00477752 0.177357` |
| `cM=2 edgeM=0 L1K3 stab=1` | 4 | `0.1875 0.1875 0.1875 416 0.000586063 0.0377799` |
| `cM=1 edgeM=0 L0K2 stab=0` | 4 | `0.375 0.375 0.375 56 0.0368586 0.662667` |

Heavier local-solve check:

| Case | Solver | Wall time | Result |
|---|---|---:|---|
| `cM=1 edgeM=3 L1K3 stab=1` | serial `MHM-RAD_poly_opt.edp` | `21.16s` | `0.375 0.046875 0.046875 896 4.70327e-06 0.00112134` |
| `cM=1 edgeM=3 L1K3 stab=1` | parallel, `np=4` | `8.85s` | `0.375 0.046875 0.046875 896 4.70327e-06 0.00112134` |
# MHMUN-RAD_Parallel
# MHMUN-RAD_Parallel
