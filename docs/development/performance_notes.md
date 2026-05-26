# FUSS Performance Optimization Notes

> Target files: `src/lib/numerics/fluxes/Mod_Fluxes.f90`, `Lib_Diffusive.f90`, and related routines

---

## Priority 1 — Cache material properties per cell instead of per face

`co_Kappa` (thermal conductivity from the material table) is called inside
`Diffusive_Flux` **once per face**, computing a face-average temperature and
performing a linear table interpolation.  For a cell shared by 6 faces, this
means the same cell-center state contributes to up to 6 separate interpolation
calls.

### Suggested fix

Precompute `kappa` at cell centres before the face loops:

```fortran
real(R8) :: kappa_cell(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc)

!$omp do collapse(3)
do k = 1-gc, n(3)+gc
do j = 1-gc, n(2)+gc
do i = 1-gc, n(1)+gc
  call co_Kappa( matID(i,j,k), T(i,j,k), kappa_cell(i,j,k) )
end do; end do; end do
```

Then average the two cell values at the face instead of interpolating from
the face-average temperature:

```fortran
kappa = 0.5d0 * ( kappa_cell(i,j,k) + kappa_cell(i+1,j,k) )
```

**Expected gain:** 10–20% (eliminates repeated table lookups in the innermost
loop).

---

## Priority 2 — Replace scalar stencil arguments with array slices

`Diffusive_Flux` currently receives the temperature stencil as ten individual
scalar arguments (`T1 … T10`).  This prevents the compiler from seeing the
full array and inhibits auto-vectorization of the face loop.

### Suggested fix

Pass the full temperature array and index coordinates:

```fortran
subroutine Diffusive_Flux ( matID, normal, area, T, M, Res, &
                             i, j, k, aa, bb, cc, n, gc )
  real(R8), intent(in)    :: T(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc)
  real(R8), intent(inout) :: Res(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc)
  integer,  intent(in)    :: i, j, k, aa, bb, cc
```

The compiler can then vectorize the `j`-`k` loop with `i` fixed, or the
calling loop can be reshaped for better SIMD utilization.

**Expected gain:** 5–15% (vectorization of inner loop).

---

## Priority 3 — Fuse direction loops

`Fluxes_blk` currently computes direction-1, direction-2, and direction-3
fluxes in separate triple loops.  Each pass streams through the full `T` and
`Res` arrays, causing redundant cache misses when the arrays do not fit in
L2/L3.

### Suggested fix

Merge all three direction sweeps into a single pass using a direction index:

```fortran
do dir = 1, ndir
  !$omp do collapse(2)
  do k = 1, n(3)
  do j = 1, n(2)
  do i = 1, n(dir) - 1
    call Diffusive_Flux( ..., dir )
  enddo; enddo; enddo
end do
```

**Expected gain:** 5–10% (improved cache reuse).

---

## Priority 4 — AoS → SoA for metric types

The derived types used for mesh metrics store small arrays inside structs:

```fortran
type :: FUSS_vector_3D_type
  real(R8) :: c(3)         ! dl(i,j,k)%c(d) — Array of Structures
end type
type :: FUSS_tensor_3D_type
  real(R8) :: c(3,3)       ! M(i,j,k)%c(d1,d2) — Array of Structures
end type
```

This layout causes stride-3 (or stride-9) access when sweeping over cells
for a single metric component, which prevents SIMD vectorization.  The
optimal layout is Structure of Arrays:

```fortran
real(R8) :: dl(3, 1-gc:n1+gc, 1-gc:n2+gc, 1-gc:n3+gc)    ! dl(d,i,j,k)
real(R8) :: M (3, 3, 1-gc:n1+gc, 1-gc:n2+gc, 1-gc:n3+gc) ! M(d1,d2,i,j,k)
```

### Impact

This is a large refactor affecting mesh setup, I/O, boundary conditions,
and all flux routines.  Best addressed incrementally, starting with the
metric tensor `M` inside `Diffusive_Flux`, where the matrix–vector product
`matmul(Gradient, M)` is the dominant operation.

**Expected gain:** 5–10% (better vectorization on the metric matmul).

---

## Priority 5 — Compiler flags

### Currently used (Release)

`-O3`, `-funroll-loops`, `-ffast-math`, `-march=native`, `-finline-functions`

### Additions to consider

| Flag (gfortran) | Effect |
|---|---|
| `-flto` | Link-time optimization — enables cross-module inlining of small routines such as `co_Kappa` and `Tangential_Gradient`. |
| `-fopt-info-vec-missed` | Diagnostic: reports what the compiler fails to vectorize — useful for identifying remaining bottlenecks. |
| `-fopenmp-simd` | Enables `!$omp simd` directives for explicit SIMD hints without full OpenMP threading overhead. |

**Expected gain:** 5–10% from `-flto` alone (cross-module inlining of
table-interpolation and gradient routines).

---

## Summary

| # | Action | Effort | Expected gain |
|---|---|---|---|
| 1 | Cache `kappa` per cell before face loops | small | 10–20% |
| 2 | Replace scalar stencil args with array + indices | medium | 5–15% |
| 3 | Fuse direction loops | small | 5–10% |
| 4 | AoS → SoA for metric types | large | 5–10% |
| 5 | Enable `-flto` | trivial | 5–10% |
