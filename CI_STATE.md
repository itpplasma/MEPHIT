# Owned CI controller handoff — 21 October 2026

[PR 75](https://github.com/itpplasma/MEPHIT/pull/75) landed through the normal
squash merge queue as a03c1293aa1d91edfc00277717f4bdeedcb6c32e at17:52:50UTC.
Authoritative merge-group run39643825414/job119374035317 passed on owned
Ubuntu24 runner7601,17:49:37–17:52:21UTC. Native dependency setup and `make`
remain unchanged, including explicit candidate LIBNEO_REF support. This
workflow builds; it does not invent extra numerical test coverage.

All PR events skip before allocation. Merge groups and direct default-branch
nonqueue pushes validate on owned Ubuntu24 capacity (two CPUs/16 GiB,
CMake parallelism two, one BLAS/OpenMP thread). The existing required
`Build and test` context and 30-minute queue deadline remain.

Existing hygiene [PR 74](https://github.com/itpplasma/MEPHIT/pull/74) was
rebased onto that actual main and sealed at
58d26f4bd4773940a66514937a6c992c5070c9d9. Its entire independent diff is an
11-line unconditional reusable-workflow caller; both shared fields pin
104634c01aa669d542e8a26a34a50b6471fc19ce. Actionlint, diff checks and configured
hooks passed. Exact-head PR CI/hygiene runs39648094189/39648094634 skipped.

The controller authorized its normal ready/squash queue after confirming
the 11-line scope and native75 main. It is now queued; merge-group runs
39648827216 and39648826010 use synthetic head
ac170f457527a9dd9c2836549fe4a6e42a5f9641. Await actual native/hygiene results;
do not equate queue admission with a merge. No source/native changes are
included. The global CI owner retains ruleset ownership. No new or replacement
runs are authorized by this handoff.
