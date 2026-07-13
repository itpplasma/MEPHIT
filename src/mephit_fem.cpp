#include "mephit_util.h"
#include "mephit_fem.h"
#ifdef USE_MFEM
#include "mfem.hpp"
#include "magnetic_differential_equation.h"
#endif  // USE_MFEM
#include <gmsh.h>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <cmath>
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

double approx_ellipse_circumf(double a, double b)
{
  const double h = (a - b) * (a - b) / ((a + b) * (a + b));
  const double fac = 1.0 + h/4 * (1.0 + h/16 * (1.0 + h/4 * (1 + 25*h/64)));
  return fac * M_PI * std::fabs(a + b);
}

// written with support from ChatGPT
extern "C" void FEM_triangulate_external_gmsh(const int npt,
                                              const double *bdry_R,
                                              const double *bdry_Z,
                                              const char *fname)
{
  const double outer_border_refinement = 0.125, outer_box_scale = 2.0;
  double R_min, R_max, R_mid, R_rad, Z_min, Z_max, Z_mid, Z_rad;

  R_min = HUGE_VAL, R_max = 0.0, Z_min = HUGE_VAL, Z_max = -HUGE_VAL;
  for (int i = 0; i < npt; i++) {
    R_min = std::min(R_min, bdry_R[i]);
    R_max = std::max(R_max, bdry_R[i]);
    Z_min = std::min(Z_min, bdry_Z[i]);
    Z_max = std::max(Z_max, bdry_Z[i]);
  }
  R_mid = 0.5 * (R_max + R_min);
  R_rad = 0.5 * (R_max - R_min) * outer_box_scale;
  Z_mid = 0.5 * (Z_max + Z_min);
  Z_rad = 0.5 * (Z_max - Z_min) * outer_box_scale;
  double size_max = approx_ellipse_circumf(R_rad, Z_rad) / (npt * outer_border_refinement);
  double dist_max = 0.5 * std::min(R_max - R_min, Z_max - Z_min) * (outer_box_scale - 1.0);

  gmsh::initialize();
  // suppress info output, but keep warnings and errors
  gmsh::option::setNumber("General.Verbosity", 3);
  gmsh::option::setNumber("Mesh.Optimize", 1);
  gmsh::option::setNumber("Mesh.OptimizeNetgen", 1);
  gmsh::option::setNumber("Mesh.MeshSizeExtendFromBoundary", 0);
  gmsh::option::setNumber("Mesh.MeshSizeFromPoints", 0);
  gmsh::option::setNumber("Mesh.MeshSizeFromCurvature", 0);
  gmsh::model::add("ext");

  // import inner boundary into gmsh and calculate mesh size auxiliaries
  double min_edge_length = HUGE_VAL, max_edge_length = 0.0;
  std::vector<double> edge_length(npt);
  for (int i = 0; i < npt; i++) {
    edge_length[i] = std::hypot(
      bdry_R[(i + 1) % npt] - bdry_R[i],
      bdry_Z[(i + 1) % npt] - bdry_Z[i]
    );
    min_edge_length = std::min(min_edge_length, edge_length[i]);
    max_edge_length = std::max(max_edge_length, edge_length[i]);
  }
  std::vector<int> point(npt);
  for (int i = 0; i < npt; i++) {
    point[i] = gmsh::model::geo::addPoint(
      bdry_R[i],
      bdry_Z[i],
      0.0
    );
  }
  std::vector<int> line(npt);
  for (int i = 0; i < npt; i++) {
    line[i] = gmsh::model::geo::addLine(
      point[i],
      point[(i + 1) % npt]
    );
  }
  int inner_loop = gmsh::model::geo::addCurveLoop(line);

  int p0 = gmsh::model::geo::addPoint(R_mid, Z_mid, 0.0);
  int p1 = gmsh::model::geo::addPoint(R_mid + R_rad, Z_mid, 0.0);
  int p2 = gmsh::model::geo::addPoint(R_mid, Z_mid + Z_rad, 0.0);
  int p3 = gmsh::model::geo::addPoint(R_mid - R_rad, Z_mid, 0.0);
  int p4 = gmsh::model::geo::addPoint(R_mid, Z_mid - Z_rad, 0.0);
  int e1 = gmsh::model::geo::addEllipseArc(p1, p0, p1, p2);
  int e2 = gmsh::model::geo::addEllipseArc(p2, p0, p1, p3);
  int e3 = gmsh::model::geo::addEllipseArc(p3, p0, p1, p4);
  int e4 = gmsh::model::geo::addEllipseArc(p4, p0, p1, p1);
  int outer_loop = gmsh::model::geo::addCurveLoop({e1, e2, e3, e4});

  // distance field yields distance from boundaries
  int dist_field = gmsh::model::mesh::field::add("Distance");
  gmsh::model::mesh::field::setNumbers(dist_field, "CurvesList", {inner_loop});
  gmsh::model::mesh::field::setNumber(dist_field, "Sampling", npt);
  // threshold field interpolates linearly from size_min to size_max
  // from distance dist_min to dist_max
  int thresh = gmsh::model::mesh::field::add("Threshold");
  gmsh::model::mesh::field::setNumber(thresh, "InField", dist_field);
  gmsh::model::mesh::field::setNumber(thresh, "SizeMin", max_edge_length);
  gmsh::model::mesh::field::setNumber(thresh, "SizeMax", size_max);
  gmsh::model::mesh::field::setNumber(thresh, "DistMin", size_max);
  gmsh::model::mesh::field::setNumber(thresh, "DistMax", dist_max);
  gmsh::model::mesh::field::setAsBackgroundMesh(thresh);

  // outer boundary first, then hole defined by inner boundary
  gmsh::model::geo::addPlaneSurface({outer_loop, inner_loop});
  gmsh::model::geo::synchronize();
  gmsh::model::mesh::generate(2);
  gmsh::model::mesh::optimize("Netgen");
  gmsh::option::setNumber("Mesh.MshFileVersion", 2.2);
  gmsh::write(fname);
  gmsh::finalize();
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
