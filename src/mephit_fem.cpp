#include "mephit_fem.h"
#ifdef USE_MFEM
#include "mfem.hpp"
#include "magnetic_differential_equation.h"
#endif  // USE_MFEM
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <cstdio>
#include <vector>
#include <map>

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

typedef bg::model::d2::point_xy<double> RZ_point;
typedef bg::model::box<RZ_point> RZ_box;
typedef std::pair<RZ_box, int> RZ_index;

static bgi::rtree< RZ_index, bgi::rstar<16> > RZ_tree;
static std::vector<int> query_results;

extern "C" void Rtree_init(int ntri, double *tri_bb)
{
  RZ_tree.clear();
  for (int k = 0; k < ntri; ++k) {
    // construct bounding box object
    RZ_box bb = RZ_box(RZ_point(tri_bb[4 * k], tri_bb[4 * k + 1]),
                       RZ_point(tri_bb[4 * k + 2], tri_bb[4 * k + 3]));
    RZ_tree.insert(std::make_pair(bb, k + 1));
  }
}

extern "C" void Rtree_query(double R, double Z, int *result_size, int **results)
{
  query_results.clear();
  for (auto it = RZ_tree.qbegin(bgi::contains(RZ_point(R, Z))); it != RZ_tree.qend(); ++it) {
    query_results.push_back(it->second);
  }
  *result_size = static_cast<int>(query_results.size());
  *results = query_results.data();
}

#ifdef USE_MFEM

typedef std::map<std::pair<double, double>, size_t> points_2D;

/*
void count_points_2D(points_2D& assoc, const double R, const double Z)
{
  auto const key = std::make_pair(R, Z);
  auto it = assoc.find(key);
  if (assoc.end() == it) {
    assoc.insert(it, std::make_pair(key, static_cast<size_t>(1)));
  } else {
    it->second += 1;
  }
}
*/

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
    points_2D points_lhs, points_rhs;
    // construct functors for real and imaginary part
    // (due to MFEM reasons, it has to be evaluated twice)
    auto const f_r = [MDE_inhom, &points_lhs](const double R, const double Z) {
      double f[2];
      MDE_inhom(R, Z, reinterpret_cast<complex_double *>(f));
      /* count_points_2D(points_lhs, R, Z); */
      return f[0];
    };
    auto const f_i = [MDE_inhom, &points_lhs](const double R, const double Z) {
      double f[2];
      MDE_inhom(R, Z, reinterpret_cast<complex_double *>(f));
      /* count_points_2D(points_lhs, R, Z); */
      return f[1];
    };
    auto const h_phi = [unit_B0, &points_rhs](const double R, const double Z) {
      double h[3];
      unit_B0(R, Z, h);
      /* count_points_2D(points_rhs, R, Z); */
      return h[2] / R;
    };
    auto const h_t = [unit_B0, &points_rhs](const double R, const double Z) {
      double h[3];
      unit_B0(R, Z, h);
      /* count_points_2D(points_rhs, R, Z); */
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
    /*
    FILE *file;
    file = fopen("lhs.txt", "w");
    if (file != nullptr) {
      for (auto const& [key, value] : points_lhs) {
        fprintf(file, "%.16e %.16e %lu\n", key.first, key.second, value);
      }
      fclose(file);
    }
    file = fopen("rhs.txt", "w");
    if (file != nullptr) {
      for (auto const& [key, value] : points_rhs) {
        fprintf(file, "%.16e %.16e %lu\n", key.first, key.second, value);
      }
      fclose(file);
    }
    */
  } catch (...) {
    return 1;
  }
  return 0;
}

#endif  // USE_MFEM
