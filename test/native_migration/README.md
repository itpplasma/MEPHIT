# Matched migration measurements

Build the parent revision and this branch independently with the same compiler,
Release flags, `BLA_VENDOR`, libneo/Triangle pins and `WITH_MFEM=OFF`. Build
`test_native_sparse` on this branch. For the original library, compile the same
`test/test_native_sparse.f90` against its generated `sparse_mod.mod` and
`libmephit.so`, preserving the same Fortran flags. Run the `matched` mode three
times in fresh processes. It performs 1000 iterations of independent known-x
real/complex scalar/two-RHS solves, with explicit factor/solve/free (`iopt=0`)
and method 3 on both libraries. The new regressions also test historical defects
that cannot be used as an equivalent successful legacy timing workload.

Run each build's registered `edge_cylinder` once, with its original scientific
tolerances and identical `OMP_NUM_THREADS=2`, `OPENBLAS_NUM_THREADS=1`. This
includes initial meshing, iteration and the separate sheet phase. Retain the
scientific check output, complete inputs and component receipts. A stopped or
failed run is not successful performance evidence.

For wall time and cumulative user/system CPU, wrap the command with
`/usr/bin/time`. Its maximum-RSS field is not an aggregate process-tree peak.
The small observer below records each actual application/FreeFem process's own
Linux post-exec `VmHWM` and self CPU at normal exit:

```sh
cc -shared -fPIC -O2 test/native_migration/measure_self.c -o measure_self.so
mkdir -p receipts
MEPHIT_METRIC_DIR="$PWD/receipts" LD_PRELOAD="$PWD/measure_self.so" \
  OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 \
  /usr/bin/time -v ctest --test-dir build -R '^edge_cylinder$' --output-on-failure -j2
```

The observer writes one JSON file per PID, including executable/parent PID.
Identify the MEPHIT computation child separately from its supervisor, and the
FreeFem daemon separately for the iteration and sheet phases. A forked MEPHIT
child inherits the near-start wall timestamp; its CPU is `RUSAGE_SELF`, and its
memory is its own `/proc/self/status` high-water. It does not use Python
fork/exec `wait4` RSS, which can inherit an unrelated parent high-water floor.
Neither sum nor maximum of individual process peaks is a simultaneous aggregate
peak; keep aggregate memory unknown unless separately measured. Ratios with a
zero baseline or below timer resolution remain unavailable.

The original default launcher may select the installed MPI/MUMPS interpreter;
the new default selects the exact configured serial interpreter and native
module. Record these executable/backend changes. The primary comparison is the
whole default workflow before and after migration, rather than isolated factor
performance. A separately provisioned minimal interpreter/SDK proves runtime
independence without conflating its build with that comparison.

A FreeFem 4.17 prefix without external sparse solver libraries was built from
commit `754dae8ae8f858784ad46b2f46dfd7862829cb5a`. Its core configure guards
include `--disable-system-umfpack --disable-umfpack --without-petsc
--without-petsc_complex --without-mpi`, with MUMPS, SuperLU and automatic
third-party downloads disabled. The exact tested provisioning recipe and
artifact hashes accompany the migration results. The interpreter may retain
dense BLAS/LAPACK, ARPACK, HDF5 and transitive MPI dependencies; the guarantee
concerns external sparse factors. Configure both the explicit executable and
matching SDK tool paths; the launcher resolves no substitute at runtime.
