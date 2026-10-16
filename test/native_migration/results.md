# Native migration validation and resource results

These are measured whole-workflow costs, not solver leadership claims. The
public cylinder, generated physical inputs and scientific checker are identical
between configurations. The baseline is original `43876ad`; native code is
`83674ec`, with FortSparse `a96033e`. GNU 16.2.1 Release builds use OpenBLAS
0.3.34, the same libneo `d30174f` and Triangle `dd24af5` artifacts, on one Ryzen
9 5950X host. Requested threads are `OMP_NUM_THREADS=2`,
`OPENBLAS_NUM_THREADS=1`; other agent CPU jobs were held during measurements.
See [machine-readable measurements](results.json) and the
[measurement recipe](README.md).

The old default uses UMFPACK in FDM and stock MPI/MUMPS FreeFem. The new default
uses native FDM factors and stock serial FreeFem with the native module. Both
builds have `WITH_MFEM=OFF`. The separately provisioned minimal interpreter uses
another FreeFem build and `WITH_MFEM=ON`; its costs have those additional changes
and are not an isolated backend comparison.

| Case | Repetitions | Legacy wall / CPU s | Native wall / CPU s | CPU change | Status |
| --- | ---: | ---: | ---: | ---: | --- |
| Matched sparse probe, median | 3 each | 0.06 / 0.064460 | 0.01 / 0.018844 | -0.045616 s (-70.77%) | both pass |
| Full public edge cylinder | 1 each | 89.22 / 92.99 | 110.42 / 112.42 | +19.43 s (+20.89%) | both pass |
| Minimal-runtime cylinder, MFEM enabled | 1 | separate configuration | 119.48 / 121.69 | not isolated | pass |

The full cylinder wall change is +21.20 s (+23.76%). Wall times include process
startup; the short wall figures have 0.01 s timer resolution. The probe performs
1000 iterations of real/complex scalar/two-RHS known-solution systems of order
3, with the unchanged 2e-14 gate. Cylinder scale is 154 surfaces, 63989 core
mesh points and a Maxwell matrix of order 194499. Its unchanged scientific
gates are 2% displacement and 3% current/field. All three configurations print
`xi/analytic=1.00017-0.00001i`, `I/analytic=1.00001-0.00020i` and field ratios
0.9808..0.9822 with maximum imaginary magnitude 0.019; the checker passes.

Peak RSS below is each process's own post-exec Linux VmHWM, in KiB. Aggregate
simultaneous memory is **unknown**; peaks are neither summed nor presented as a
process-tree peak. Whole-workflow CPU above is cumulative user+system from
`/usr/bin/time`; component CPU below is `RUSAGE_SELF`.

| Process | Legacy CPU s | Native CPU s | CPU change | Legacy peak KiB | Native peak KiB | Peak change |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Matched probe, median | 0.064460 | 0.018844 | -0.045616 (-70.77%) | 18964 | 6536 | -12428 (-65.53%) |
| Iteration MEPHIT compute child | 8.830960 | 8.748936 | -0.082024 (-0.93%) | 933380 | 930188 | -3192 (-0.34%) |
| Iteration FreeFem daemon | 6.441397 | 14.474414 | +8.033017 (+124.71%) | 300596 | 496316 | +195720 (+65.11%) |
| Sheet MEPHIT compute child | 67.741054 | 72.063516 | +4.322462 (+6.38%) | 2529548 | 2460512 | -69036 (-2.73%) |
| Sheet FreeFem daemon | 7.827760 | 16.190922 | +8.363162 (+106.84%) | 301792 | 496420 | +194628 (+64.49%) |

The minimal-runtime configuration's iteration MEPHIT/FreeFem peaks are
1273260/495832 KiB; sheet peaks are 2460884/495820 KiB. Its corresponding CPU
values are 16.194947/14.894742 s and 73.110462/16.495311 s. Its checked closure
excludes UMFPACK, CHOLMOD, SuperLU and MUMPS for libmephit, the interpreter,
BinaryIO and the shared native module. The stock interpreter retains unused
vendor DSOs; only its selected sparse execution is native. See the
[tested minimal-runtime provisioning recipe](minimal-freefem.md).

Focused stock-default tests pass 10/10; the minimal-runtime optional-MFEM build
passes 11/11. They cover known solutions, real/complex scalar/block and retained
factors, singular/error behavior, linkage, FreeFem transpose/reuse, delayed FIFO
framing with 100 actual BinaryIO cycles, MFEM CSR and the actual penalty routine,
and the existing edge refinement unit. A full physical MFEM case is not claimed.

The first native cylinder attempt (`d9d45d2`) failed after 98.49 s with FIFO EOF
before scientific checking; it has no successful performance interpretation.
An independent delayed consecutive-writer reproducer fails against the original
library with the same 0/8-byte EOF. Two vector-consumption acknowledgments fix
the framing race without changing equations, Richardson corrections or gates.
Raw failed/successful receipts, artifact hashes, full check output and runtime
closure are retained in the task's `consumer-migrations/mephit` evidence bundle;
no generated HDF5 output or private fixture is committed.

The final dependency pin is `19f9b4734a8e67faf05dac41bb19a2f8821eaaae`.
The CPU/RSS tables above remain the exact `83674ec` / `a96033e` measured
checkpoint. Later upstream differences are an explicit symmetry-normalization
helper that MEPHIT does not call, documentation, and FreeFem SDK configuration
header portability (`FreeFEM-config.h` in Ubuntu 4.13, `config.h` in 4.17).
Existing numerical factor/solve routines are unchanged. Affected consumer
checks and linkage were repeated after the pin; no new head timing is inferred.

Ubuntu CI additionally exposed packaged BinaryIO/idp files outside the compiled
FreeFem search paths. Matching-prefix resource discovery now configures those
paths explicitly for launchers and tests. Stock/minimal affected FreeFem checks
pass again; the actual official Ubuntu 4.13 BinaryIO runtime also passes the
100-cycle exchange check with its packaged idp/plugin paths. No equations or
numerical options change, and the earlier timing checkpoint remains explicit.
