#include "mephit_fem.h"
#include "mfem.hpp"
#include "magnetic_differential_equation.h"
#include <cstdio>

extern "C" int FEM_test(const char *mesh_file,
                        const int tor_mode,
                        const int n_dof,
                        complex_double *dof,
                        real_vector_field *unit_B0,
                        complex_scalar_field *MDE_inhom)
{
  int const fe_space_order = 1;
  int const quadrature_order = 2;
  try {
    // construct functors for real and imaginary part
    // (due to MFEM reasons, it has to be evaluated twice)
    auto const f_r = [MDE_inhom](const double R, const double Z) {
      double f[2];
      MDE_inhom(R, Z, reinterpret_cast<complex_double *>(f));
      return f[0];
    };
    auto const f_i = [MDE_inhom](const double R, const double Z) {
      double f[2];
      MDE_inhom(R, Z, reinterpret_cast<complex_double *>(f));
      return f[1];
    };
    auto const h_phi = [unit_B0](const double R, const double Z) {
      double h[3];
      unit_B0(R, Z, h);
      return h[2];
    };
    auto const h_t = [unit_B0](const double R, const double Z){
      double h[3];
      unit_B0(R, Z, h);
      mfem::Vector h_vec(2);
      h_vec(0) = h[0];
      h_vec(1) = h[1];
      return h_vec;
    };
    // read mesh
    mfem::Mesh mesh(mesh_file);
    // construct matrix and rhs
    // get system and rhs
    const int dim = mesh.Dimension();
    mfem::H1_FECollection u_coll(fe_space_order, dim);
    mfem::FiniteElementSpace U(&mesh, &u_coll);
    mfem::SparseMatrix Sys = mde::construct_mde_matrix(U, h_t, h_phi, tor_mode, quadrature_order);
    mfem::Vector rhs = mde::construct_mde_rhs(U, h_t, h_phi, f_r, f_i, tor_mode, quadrature_order);
    // solve system
    // construct solver
    // If MFEM was compiled with SuiteSparse, use UMFPACK to solve the system.
    mfem::UMFPackSolver umf_solver;
    umf_solver.Control[UMFPACK_ORDERING] = UMFPACK_ORDERING_METIS;
    umf_solver.SetOperator(Sys);
    mfem::Vector sol(rhs.Size());
    umf_solver.Mult(rhs, sol);
    if (n_dof != U.GetVSize()) {
      const int dof_mesh = U.GetVSize();
      fprintf(stderr, "Mismatch of vector sizes: n_dof = %i, GetVSize() -> %i.\n", n_dof, dof_mesh);
      return 1;
    }
    for (int i = 0; i < n_dof; ++i){
      reinterpret_cast<double *>(dof)[2 * i] = sol(i);
      reinterpret_cast<double *>(dof)[2 * i + 1] = sol(n_dof + i);
    }
  } catch (...) {
    return 1;
  }
  return 0;
}
