#include "mephit_native_mfem.h"
#include "fortsparse_native.h"
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#if defined(MFEM_USE_SUITESPARSE) || defined(MFEM_USE_SUPERLU) || defined(MFEM_USE_SUPERLU5) || defined(MFEM_USE_MUMPS)
#error "MEPHIT requires MFEM without external sparse factorization libraries"
#endif
#ifdef MFEM_USE_MPI
#include <mpi.h>
#endif

namespace {
struct NativeFactor {
  void *handle = nullptr;
  ~NativeFactor() { fs_native_free(&handle); }
  void check(int32_t status) const {
    if (status == 0) { return; }
    char detail[512] = {};
    fs_native_last_error(handle, detail, sizeof(detail));
    throw std::runtime_error("FortSparse native failure " + std::to_string(status) + ": " + detail);
  }
};
}

void mephit::solve_native(const mfem::SparseMatrix &matrix,
                          const mfem::Vector &rhs, mfem::Vector &solution)
{
#ifdef MFEM_USE_MPI
  int initialized = 0, size = 1;
  MPI_Initialized(&initialized);
  if (initialized) { MPI_Comm_size(MPI_COMM_WORLD, &size); }
  if (size != 1) { throw std::runtime_error("FortSparse MEPHIT adapter requires one MPI process"); }
#endif
  static_assert(sizeof(int) == sizeof(int32_t), "MFEM CSR requires 32-bit indices");
  const int n = matrix.Height();
  if (!matrix.Finalized() || n <= 0 || n == std::numeric_limits<int32_t>::max() ||
      matrix.Width() != n || rhs.Size() != n) {
    throw std::runtime_error("FortSparse MEPHIT adapter requires finalized square CSR and matching RHS");
  }
  const int *offset = matrix.GetI();
  const int *column = matrix.GetJ();
  const double *values = matrix.GetData();
  const int nnz = offset[n];
  if (offset[0] != 0 || nnz < 0 || nnz == std::numeric_limits<int32_t>::max()) { throw std::runtime_error("Invalid MFEM CSR offsets"); }
  std::vector<int32_t> pointer(n + 1, 0), row(nnz), next(n);
  std::vector<double> data(nnz);
  for (int i = 0; i < n; ++i) {
    if (offset[i] < 0 || offset[i] > offset[i + 1] || offset[i + 1] > nnz) {
      throw std::runtime_error("Invalid MFEM CSR offsets");
    }
    for (int k = offset[i]; k < offset[i + 1]; ++k) {
      if (column[k] < 0 || column[k] >= n) { throw std::runtime_error("Invalid MFEM CSR column"); }
      ++pointer[column[k] + 1];
    }
  }
  for (int j = 0; j < n; ++j) { pointer[j + 1] += pointer[j]; next[j] = pointer[j]; }
  // Ascending original rows produce sorted CSC row indices without changing values.
  for (int i = 0; i < n; ++i) {
    for (int k = offset[i]; k < offset[i + 1]; ++k) {
      const int target = next[column[k]]++;
      row[target] = i;
      data[target] = values[k];
    }
  }
  NativeFactor factor;
  factor.check(fs_native_create(0, &factor.handle));
  factor.check(fs_native_set_ordering(factor.handle, FS_NATIVE_ORDER_AMD));
  // The former MFEM UMFPACK solver enabled iterative refinement by default.
  factor.check(fs_native_set_refine(factor.handle, 1));
  factor.check(fs_native_factor_d(factor.handle, n, nnz, pointer.data(), row.data(), data.data()));
  mfem::Vector candidate(n);
  factor.check(fs_native_solve_d(factor.handle, 0, n, 1, rhs.GetData(), candidate.GetData()));
  solution = candidate;
}
