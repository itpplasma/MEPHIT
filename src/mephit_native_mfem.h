#ifndef MEPHIT_NATIVE_MFEM_H
#define MEPHIT_NATIVE_MFEM_H
#include "mfem.hpp"
namespace mephit {
// Finalized serial CSR is converted to CSC; FortSparse owns factorization.
void solve_native(const mfem::SparseMatrix &matrix, const mfem::Vector &rhs,
                  mfem::Vector &solution);
}
#endif
