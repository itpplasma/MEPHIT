# Contributing to MEPHIT

Propose changes as pull requests against `main` and run `make test` first.

## Continuous integration and merging

- Every pull request, including each layer of a stacked pull request, runs CI.
- Review a stack layer by layer. When it is ready, retarget the top layer to
  `main` so that the stack is merged as one pull request; the lower pull
  requests then show as merged.
- Merge with "Merge when ready" (or `gh pr merge --merge --auto`). The pull
  request enters the merge queue, CI runs on the exact merge result, and on
  success it lands on `main` as a merge commit. `main` is not re-tested on
  push.
- A branch does not need to be up to date with `main`. Merge `main` into a
  pull request only when it conflicts.
